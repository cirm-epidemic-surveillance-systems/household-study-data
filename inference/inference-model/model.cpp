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

// #define COMPLETE_MODEL_DRAFT

#ifdef COMPLETE_MODEL_DRAFT
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
	double date, inst_val, cumul_val;

	std::ifstream input(pathtofile.c_str());
	if (input.is_open())
	{
		// load instantaneous value of the epi curve
		for (std::string line; getline(input, line, ',');)
		{
			std::istringstream iss(line);

			// Store information in variables
			iss >> date >> inst_val;
			Date_EPI_Curve.push_back(date);
			Inst_EPI_Curve.push_back(inst_val);
		}

		// store the cumulative epi curve
		Cumul_EPI_Curve[0] = Inst_EPI_Curve[0];
		for (size_t i = 1; i < Inst_EPI_Curve.size(); ++i)
			Cumul_EPI_Curve[i] = Inst_EPI_Curve[i] + Cumul_EPI_Curve[i - 1];
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
	else // calculate the interpolated y value
	{
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
}

/* Interpolated value of the instantaneous epi curve. */
double epi_curve_inst(double t) { return interpolate_by_date(t, Inst_EPI_Curve); }

/* Interpolated value of the cumulative epi curve. */
double epi_curve_cumul(double t) { return interpolate_by_date(t, Cumul_EPI_Curve); }
#endif

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/* Assign the contact patterns between a pair of individuals. */
void Model::assign_contact_rate_to_household(Household &house, indiv from, indiv to)
{
#ifdef COMPLETE_MODEL_DRAFT
	const bool from_child = house.data(from, "age") <= 15;
	const bool from_senior = house.data(from, "age") >= 65;
	const bool from_adult = (!from_child) && (!from_senior);

	const bool to_child = house.data(from, "age") <= 15;
	const bool to_senior = house.data(from, "age") >= 65;
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
		emplace_contact_rate("ctct_adu2child");
	else if (from_adult && to_adult)
		emplace_contact_rate(REFERENCE);
	else if (from_adult && to_senior)
		emplace_contact_rate("ctct_adu2senior");

	// from is an adult
	else if (from_child && to_child)
		emplace_contact_rate("ctct_senior2child");
	else if (from_child && to_adult)
		emplace_contact_rate("ctct_senior2adult");
	else if (from_child && to_senior)
		emplace_contact_rate("ctct_senior2senior");

	else
		emplace_contact_rate(MISSING_RATE);
#else
	UNUSED(house);
	UNUSED(from);
	UNUSED(to);
	emplace_contact_rate(REFERENCE);
#endif
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/* Assign all group rates to each individual. */
void Model::assign_group_rate_to_individual(Household &house, indiv i)
{
	rateID r;
	UNUSED(house);
	UNUSED(i);

	//----- infectivity -----//
	r = "infectivity";
	emplace_group_rate(r, REFERENCE);

	//----- susceptibility -----//
	r = "susceptibility";
	emplace_group_rate(r, REFERENCE);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Model : initialization  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize infectivity profile. */
void Model::init_infec_profile()
{
	// default infectivity profile: gamma
	// otherwise user defined
	assign_template_infec_profile({{"name", "reparam_gamma"},
								   {"arg_name", "mean/var"},
								   {"arg_value", "2/1"},
								   {"arg_at", "S/S"},
								   {"discrete", "TRUE"}});
}

/* Initialize parameters and augmented variables. */
void Model::init_model_for_inference()
{
	// integer for index of introduction case
	allocate_extra_house_integer_vector("intro_case");

	// matrices for the instantaneous and cumulative infectivity profiles
	allocate_extra_house_double_matrix("inst");
	allocate_extra_house_double_matrix("cumul");

#ifdef COMPLETE_MODEL_DRAFT
	// read epi curve from csv file
	epi_curve_read_csv("epi-curve.csv");
#endif
}

/* Initialize parameters and augmented variables. */
void Model::init_model_for_simulation()
{
	init_model_for_inference();
}

/* Transmission rate. */
double Model::transmission_rate(Household &house) const
{
	return param("beta") * std::pow(param("size_ref") / house.size(), param("gamma"));
}

/* Cumulative Risk of Infection for susceptible individual i. */
double Model::loglike_healthy(Household &house) const
{
	const indiv i0 = house.extra_integer_vector("intro_case")[0];
	const double t0 = house.augm(i0, "first_pos_date") - house.augm(i0, "delay_period");
	const double tf = t0 + param("study_period");

	const double a = param("alpha");
	const double b = transmission_rate(house);

	// load pre-compute inefctivity profiles
	const auto &cumul_mtx = house.extra_double_matrix("cumul");

	double output = 0.;

	for (indiv i : house.healthy_status_list())
	{
		double roi_cumul_i = 0;
		for (indiv k : house.infected_status_list())
			roi_cumul_i += cumul_mtx[k][i] * house.contact_rate(k, i) * house.group_rate(k, "infectivity");
		roi_cumul_i *= house.group_rate(i, "susceptibility") * b;
		output -= roi_cumul_i;
	}

#ifdef COMPLETE_MODEL_DRAFT
	output -= a * (epi_curve_cumul(tf) - epi_curve_cumul(t0)) * house.n_healthy_status();
#else
	output -= a * (tf - t0) * house.n_healthy_status();
#endif

	return output;
}

/* Cumulative Risk of Infection for infectious individual i. */
double Model::loglike_infectious(Household &house) const
{
	double output = 0.;
	const indiv i0 = house.extra_integer_vector("intro_case")[0];
	const double t0 = house.augm(i0, "first_pos_date") - house.augm(i0, "delay_period");

	const auto &inst_mtx = house.extra_double_matrix("inst");
	const auto &cumul_mtx = house.extra_double_matrix("cumul");

	const double a = param("alpha");
	const double b = transmission_rate(house);

	for (indiv i : house.infected_status_list())
		if (i != i0)
		{
			// Risk of infection from other household members `k`
			double b_i = b * house.group_rate(i, "susceptibility");
			double roi_inst_i = 0.;
			double roi_cumul_i = 0.;
			for (indiv k : house.infected_status_list())
				if (cumul_mtx[k][i] != 0)
				{
					double mu_ki = house.contact_rate(k, i) * house.group_rate(k, "infectivity");
					roi_inst_i += inst_mtx[k][i] * mu_ki;
					roi_cumul_i += cumul_mtx[k][i] * mu_ki;
				}
			roi_inst_i *= b_i;
			roi_cumul_i *= b_i;

			// Risk of infection from community
			double t_infec_i = house.augm(i, "first_pos_date") - house.augm(i, "delay_period");

#ifdef COMPLETE_MODEL_DRAFT
			roi_inst_i += a * epi_curve_inst(t_infec_i);
			roi_cumul_i += a * (epi_curve_cumul(t_infec_i) - epi_curve_cumul(t0));
#else
			roi_inst_i += a;
			roi_cumul_i += a * (t_infec_i - t0 - 1);
#endif
			output += std::log(1. - std::exp(-roi_inst_i)) - roi_cumul_i;
		}

	return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store cumulative and instantaneous infectivity profiles according to an infector.
   WARNING: unpredictable behavior if i is not an infector.
   WARNING: this function should not call any mutable parameter. */
void Model::update_extra_dependecies(Household &house, indiv i)
{
	const double onset_i = house.augm(i, "first_pos_date");
	const double incub_i = house.augm(i, "delay_period");
	const double t_infec_i = onset_i - incub_i;

	auto &inst_mtx = house.extra_double_matrix("inst");
	auto &cumul_mtx = house.extra_double_matrix("cumul");

	// init new first case to i
	double new_t0 = t_infec_i;
	indiv i0 = i;

	//==============================//
	//===== Infectious Members =====//
	//==============================//

	// note: setting lut to 0 this is necessary because _houseROI_ is brute
	for (indiv k : house.infected_status_list())
		if (k != i)
		{
			inst_mtx[i][k] = 0.0;
			cumul_mtx[i][k] = 0.0;
			inst_mtx[k][i] = 0.0;
			cumul_mtx[k][i] = 0.0;

			double incub_k = house.augm(k, "delay_period");
			double t_infec_k = house.augm(k, "first_pos_date") - incub_k;
			double gentime = t_infec_k - t_infec_i;

			if (gentime > 0)
			{ //----- `i` : infector, `k` : infectee -----//
				inst_mtx[i][k] = infec_profile_inst(gentime);
				cumul_mtx[i][k] = infec_profile_cumul(gentime -1);
			}
			else if (gentime < 0)
			{ //----- `k` : infector, `i` : infectee -----//
				inst_mtx[k][i] = infec_profile_inst(-gentime);
				cumul_mtx[k][i] = infec_profile_cumul(-gentime - 1);

				//----- update first case -----//
				if (t_infec_k < new_t0)
				{
					new_t0 = t_infec_k;
					i0 = k;
				}
			}
			// case gentime == 0 is skipped
		}

	// store intro case
	house.extra_integer_vector("intro_case")[0] = i0;

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	double exposure = param("study_period") + house.augm(i0, "first_pos_date") - t_infec_i;
	for (indiv k : house.healthy_status_list())
	{
		cumul_mtx[i][k] = infec_profile_cumul(exposure);
		cumul_mtx[k][i] = 0.0;
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_house_obs_loglike(Household &house) const
{
	return loglike_healthy(house) + loglike_infectious(house);
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

/* Instantaneous Risk of Infection for individual i (for simulation). */
double Model::roi_inst_simulation(const std::vector<indiv> &infectors,
								  Household &house, indiv i, double t_infec_i) const
{
	double hh_roi_inst = 0.;

	for (indiv k : infectors)
	{
		// `k` : infector, `i` : infectee
		double incub_k = house.augm(k, "delay_period");
		double gentime = t_infec_i - house.augm(k, "first_pos_date") + incub_k;
		double mu_ki = house.contact_rate(k, i) * house.group_rate(k, "infectivity");
		hh_roi_inst += infec_profile_inst(gentime, {incub_k}) * mu_ki;
	}
	hh_roi_inst *= house.group_rate(i, "susceptibility") * transmission_rate(house);

#ifdef COMPLETE_MODEL_DRAFT
	return param("alpha") * epi_curve_inst(t_infec_i) + hh_roi_inst;
#else
	return param("alpha") + hh_roi_inst;
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Model : simulate single house  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson point process to simulate outbreak in a single household. */
void Model::simulate_single_house(Household &house)
{
	//=================================================//
	//===== GENERATION 0 : by single introduction =====//
	//=================================================//

	// pick single introduction
	indiv i0;
	if (house.n_infec_status() > 0)
		// amongst pre-existing infected in data
		i0 = house.infected_status_list().at(uniform() * house.n_infec_status());
	else
		// amongst all members
		i0 = uniform() * house.size();

	// set all to healthy except from introduction
	// assign augmented variables to all
	for (indiv i = 0; i < house.size(); ++i)
	{
		house.change_status(i, (i == i0) ? INFECTED_STATUS : HEALTHY_STATUS);
		assign_all_augm(house, i);
	}

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	// first infection occuring at time t = 0
	const double t0 = house.augm(i0, "first_pos_date") - house.augm(i0, "delay_period");
	const double tf = house.augm(i0, "first_pos_date") + param("study_period");

	// infection of secondary cases following the first case
	for (double t = t0 + 1; (t < tf) && (house.n_healthy_status() > 0); t += 1)
	{
		// store a copy of the vector of infected individuals
		const std::vector<indiv> infectors = house.infected_status_list();

		for (indiv i : house.healthy_status_list())
			if (uniform() >= exp(-roi_inst_simulation(infectors, house, i, t)))
			{
				house.set_augm_value(i, "first_pos_date", t + house.augm(i, "delay_period"));
				house.change_status(i, INFECTED_STATUS);
			}
	}

	// format infected with unseen onset as healthy
	for (indiv i : house.infected_status_list())
		if (house.augm(i, "first_pos_date") >= tf)
			house.change_status(i, HEALTHY_STATUS);

	// format all augmented variables of susceptibles
	for (indiv i : house.healthy_status_list())
		clear_all_augm(house, i);
}
