/*
 * model.cpp
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////  * *** * AUGMENTED ONSET - SINGLE INTRODUCTION * *** *  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

// use the EPI curve for the community contribution
std::vector<double> Date_EPI_Curve = {};
std::vector<double> Inst_EPI_Curve = {};
std::vector<double> Cumul_EPI_Curve = {};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////  epi curve  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read EPI curve from csv file. */
void epi_curve_read_csv(std::string pathtofile)
{
	std::ifstream input(pathtofile.c_str());
	if (input.is_open())
	{
		// load instantaneous value of the epi curve
		std::string line;
		while (getline(input, line))
		{
			std::stringstream ss(line);
			std::string item;

			// Store information in variables
			if (std::getline(ss, item, ','))
				Date_EPI_Curve.push_back(std::stod(item));
			if (std::getline(ss, item, ','))
				Inst_EPI_Curve.push_back(std::stod(item));
		}

		// store the cumulative epi curve
		Cumul_EPI_Curve.push_back(Inst_EPI_Curve.front());
		for (size_t i = 1; i < Inst_EPI_Curve.size(); ++i)
			Cumul_EPI_Curve.push_back(Inst_EPI_Curve[i] + Cumul_EPI_Curve[i - 1]);
	}
	else
		tools::error("Impossible to read file: " + pathtofile + ".");
}

/* Linear interpolation by date. */
double interpolate_by_date(double x0, const std::vector<double> &y)
{
	const std::vector<double> &x = Date_EPI_Curve; // alias

	if (x0 <= x.front()) // saturate to lower value
		return y.front();
	else if (x0 >= x.back()) // saturate to upper value
		return y.back();

	// else calculate the interpolated y value
	// bisection
	size_t i = 0;			 // lower bound
	size_t u = x.size() - 1; // upper bound
	while (u - i > 1)
	{
		size_t mid = i + (u - i) / 2;
		if (x[mid] <= x0)
			i = mid;
		else
			u = mid;
	}
	return y[i] + (x0 - x[i]) * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
}

/* Interpolated value of the instantaneous epi curve. */
double epi_curve_inst(double t) { return interpolate_by_date(t, Inst_EPI_Curve); }

/* Interpolated value of the cumulative epi curve. */
double epi_curve_cumul(double t) { return interpolate_by_date(t, Cumul_EPI_Curve); }

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/* Assign the contact patterns between a pair of individuals. */
void Model::assign_contact_rate_to_household(Household &house, indiv from, indiv to)
{
	Individual &member_from = house.get_member(from);
	const bool from_child = member_from.data("age") <= 15;
	const bool from_senior = member_from.data("age") >= 65;
	const bool from_adult = (!from_child) && (!from_senior);

	Individual &member_to = house.get_member(to);
	const bool to_child = member_to.data("age") <= 15;
	const bool to_senior = member_to.data("age") >= 65;
	const bool to_adult = (!to_child) && (!to_senior);

	// from is a child
	if (from_child && to_child)
		emplace_contact_rate("ctct_child2child");
	else if (from_child && to_adult)
		emplace_contact_rate("ctct_child2adult");
	else if (from_child && to_senior)
		emplace_contact_rate("ctct_child2senior");

	// from is an adult
	else if (from_adult && to_child)
		emplace_contact_rate("ctct_adult2child");
	else if (from_adult && to_adult)
		emplace_contact_rate(REFERENCE);
	else if (from_adult && to_senior)
		emplace_contact_rate("ctct_adult2senior");

	// from is an adult
	else if (from_senior && to_child)
		emplace_contact_rate("ctct_senior2child");
	else if (from_senior && to_adult)
		emplace_contact_rate("ctct_senior2adult");
	else if (from_senior && to_senior)
		emplace_contact_rate("ctct_senior2senior");

	else
		emplace_contact_rate(MISSING_RATE);
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/* Assign all group rates to each individual. */
void Model::assign_group_rate_to_individual(Individual &member)
{
	rateID r;
	const bool is_child = member.data("age") <= 15;
	const bool is_senior = member.data("age") >= 65;
	const bool is_adult = (!is_child) && (!is_senior);

	//----- community -----//
	r = "community";
	if (is_child)
		emplace_group_rate(r, "com_child");
	else if (is_adult)
		emplace_group_rate(r, REFERENCE);
	else if (is_senior)
		emplace_group_rate(r, "com_senior");
	else
		emplace_group_rate(r, MISSING_RATE);

	//----- infectivity -----//
	r = "infectivity";
	if (is_child)
		emplace_group_rate(r, "inf_child");
	else if (is_adult)
		emplace_group_rate(r, REFERENCE);
	else if (is_senior)
		emplace_group_rate(r, "inf_senior");
	else
		emplace_group_rate(r, MISSING_RATE);

	//----- susceptibility -----//
	r = "susceptibility";
	if (is_child)
		emplace_group_rate(r, "sus_child");
	else if (is_adult)
		emplace_group_rate(r, REFERENCE);
	else if (is_senior)
		emplace_group_rate(r, "sus_senior");
	else
		emplace_group_rate(r, MISSING_RATE);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Model : initialization  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize parameters and augmented variables. */
void Model::init_model_for_inference()
{
	tools::error("Attempt to call inference method of pure simulation model.");
}

/* Initialize parameters and augmented variables. */
void Model::init_model_for_simulation()
{
	// read epi curve from csv file
	epi_curve_read_csv("simulation/epi-curve.csv");
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Model : risk of infection for simulation  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Transmission rate. */
double Model::transmission_rate(Household &house) const
{
	return param("beta") * std::pow(param("size_ref") / house.size(), param("gamma"));
}

/* Instantaneous Risk of Infection for individual i (for simulation). */
double Model::roi_inst_simulation(const std::vector<indiv> &infectors,
								  Household &house,
								  Individual &member_i,
								  double t_infec_i) const
{
	indiv i = member_i.hhindex();
	double hh_roi_inst = 0.;

	const double a_i = param("alpha") * member_i.group_rate("community");

	// `k` : infector, `i` : infectee
	for (indiv k : infectors)
	{
		Individual &member_k = house.get_member(k);
		double mu_ki = house.contact_rate(k, i) * member_k.group_rate("infectivity");
		double gentime = t_infec_i - member_k.augm("infec_date");
		hh_roi_inst += mu_ki * infec_profile_inst(member_k, gentime);
	}
	hh_roi_inst *= member_i.group_rate("susceptibility") * transmission_rate(house);

	return a_i * epi_curve_inst(t_infec_i) + hh_roi_inst;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Model : simulate single house  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson point process to simulate outbreak in a single household. */
void Model::simulate_single_house(Household &house)
{
	const double Dt = param("delta_t");
	const double t0 = Date_EPI_Curve.front();
	const double tf = Date_EPI_Curve.back();

	// pre-process household members
	for (Individual &member : house)
	{
		house.change_status(member, HEALTHY_STATUS); // set member to healthy
		assign_all_augm(member);					 // assign augmented variables
	}

	// infection household members
	for (double t = t0; (t < tf) && (house.n_healthy_status() > 0); t += Dt)
	{
		// store a copy of the vector of infected individuals
		const std::vector<indiv> infectors = house.get_infected_status_hhindexes();

		for (Individual &member : house.healthy())
			if (uniform() >= std::exp(-roi_inst_simulation(infectors, house, member, t) * Dt))
			{
				member.set_augm_value("infec_date", t);
				house.change_status(member, INFECTED_STATUS);

				// assign recovery date
				// choosing time at which cumulative infectivity profile ~ 1
				double s = member.augm("ip_mean");
				for (; s < 50; s += 2 * Dt)
					if (1.e-5 > 1. - infec_profile_cumul(member, s))
						break;
				member.set_augm_value("recov_date", t + s);
			}
	}

	// format all augmented variables of susceptibles
	for (Individual &member : house.healthy())
		clear_all_augm(member);
}
