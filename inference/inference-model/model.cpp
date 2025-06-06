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

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/* Assign the contact patterns between a pair of individuals.
 * [PLACEHOLDER]
 */
void Model::assign_contact_rate_to_household(Household &, indiv, indiv)
{
	emplace_contact_rate(REFERENCE);
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/* Assign all group rates to each individual.
 * [PLACEHOLDER]
 */
void Model::assign_group_rate_to_individual(Individual &)
{
	rateID r;

	//----- community -----//
	r = "community";
	emplace_group_rate(r, REFERENCE);

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
	assign_cache_house_integer_vector("intro_case");

	// matrices for the instantaneous and cumulative infectivity profiles
	assign_cache_house_double_matrix("inst");
	assign_cache_house_double_matrix("cumul");
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
	const indiv i0 = house.cache_integer_vector("intro_case")[0];
	const Individual &member_0 = house.get_member(i0);
	const double t0 = member_0.augm("first_pos_date") - member_0.augm("delay_period");
	const double tf = t0 + param("study_period");

	const double a = param("alpha");
	const double b = transmission_rate(house);

	// load pre-compute inefctivity profiles
	const auto &cumul_mtx = house.cache_double_matrix("cumul");

	double output = 0.;

	for (Individual &member_i : house.healthy())
	{
		indiv i = member_i.hhindex();
		double roi_cumul_i = 0;
		for (Individual &member_k : house.infected())
		{
			indiv k = member_k.hhindex();
			roi_cumul_i += cumul_mtx[k][i] * house.contact_rate(k, i) * member_k.group_rate("infectivity");
		}
		roi_cumul_i *= member_i.group_rate("susceptibility") * b;
		output -= roi_cumul_i;
		output -= a * (tf - t0) * member_i.group_rate("community");
	}

	return output;
}

/* Cumulative Risk of Infection for infectious individual i. */
double Model::loglike_infectious(Household &house) const
{
	double output = 0.;
	const indiv i0 = house.cache_integer_vector("intro_case")[0];
	Individual &member_0 = house.get_member(i0);
	const double t0 = member_0.augm("first_pos_date") - member_0.augm("delay_period");

	const auto &inst_mtx = house.cache_double_matrix("inst");
	const auto &cumul_mtx = house.cache_double_matrix("cumul");

	const double a = param("alpha");
	const double b = transmission_rate(house);

	for (Individual &member_i : house.infected())
	{
		indiv i = member_i.hhindex();
		double b_i = b * member_i.group_rate("susceptibility");

		if (member_i != member_0)
		{
			// Risk of infection from other household members `k`
			double roi_inst_i = 0.;
			double roi_cumul_i = 0.;
			for (Individual & member_k : house.infected())
			{
				indiv k = member_k.hhindex();
				if (cumul_mtx[k][i] != 0)
				{
					double mu_ki = house.contact_rate(k, i) * member_k.group_rate("infectivity");
					roi_inst_i += inst_mtx[k][i] * mu_ki;
					roi_cumul_i += cumul_mtx[k][i] * mu_ki;
				}
			}
			roi_inst_i *= b_i;
			roi_cumul_i *= b_i;

			// Risk of infection from community
			double t_infec_i = member_i.augm("first_pos_date") - member_i.augm("delay_period");
			roi_inst_i += a * member_i.group_rate("community");
			roi_cumul_i += a * (t_infec_i - t0 - 1) * member_i.group_rate("community");
			output += std::log(1. - std::exp(-roi_inst_i)) - roi_cumul_i;
		}
	}
	return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store cumulative and instantaneous infectivity profiles according to an infector.
   WARNING: unpredictable behavior if i is not an infector.
   WARNING: this function should not call any mutable parameter. */
void Model::update_cache_dependecies(Household &house, indiv i)
{
	Individual &member_i = house.get_member(i);
	const double onset_i = member_i.augm("first_pos_date");
	const double incub_i = member_i.augm("delay_period");
	const double t_infec_i = onset_i - incub_i;

	auto &inst_mtx = house.cache_double_matrix("inst");
	auto &cumul_mtx = house.cache_double_matrix("cumul");

	// init new first case to i
	double new_t0 = t_infec_i;
	indiv i0 = i;

	//==============================//
	//===== Infectious Members =====//
	//==============================//

	// note: setting lut to 0 this is necessary because _houseROI_ is brute
	for (Individual & member_k : house.infected())
		if (member_k != member_i)
		{
			indiv k = member_k.hhindex();
			inst_mtx[i][k] = 0.0;
			cumul_mtx[i][k] = 0.0;
			inst_mtx[k][i] = 0.0;
			cumul_mtx[k][i] = 0.0;

			double incub_k = member_k.augm("delay_period");
			double t_infec_k = member_k.augm("first_pos_date") - incub_k;
			double gentime = t_infec_k - t_infec_i;

			if (gentime > 0)
			{ //----- `i` : infector, `k` : infectee -----//
				inst_mtx[i][k] = infec_profile_inst(member_i, gentime);
				cumul_mtx[i][k] = infec_profile_cumul(member_i, gentime - 1);
			}
			else if (gentime < 0)
			{ //----- `k` : infector, `i` : infectee -----//
				inst_mtx[k][i] = infec_profile_inst(member_k, -gentime);
				cumul_mtx[k][i] = infec_profile_cumul(member_k, -gentime - 1);

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
	house.cache_integer_vector("intro_case")[0] = i0;
	Individual &member_0 = house.get_member(i0);

	//===============================//
	//===== Susceptible Members =====//
	//===============================//

	const double exposure = param("study_period") + member_0.augm("first_pos_date") - t_infec_i;
	const double cumul_ip = infec_profile_cumul(member_i, exposure);
	for (indiv k : house.get_healthy_status_hhindexes())
	{
		cumul_mtx[i][k] = cumul_ip;
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
								  Household &house, 
								  Individual &member_i, 
								  double t_infec_i) const
{
	indiv i = member_i.hhindex();
	double hh_roi_inst = 0.;

	for (indiv k : infectors)
	{
		Individual &member_k = house.get_member(k);
		// `k` : infector, `i` : infectee
		double incub_k = member_k.augm("delay_period");
		double gentime = t_infec_i - member_k.augm("first_pos_date") + incub_k;
		double mu_ki = house.contact_rate(k, i) * member_k.group_rate("infectivity");
		hh_roi_inst += infec_profile_inst(member_k, gentime) * mu_ki;
	}
	hh_roi_inst *= member_i.group_rate("susceptibility") * transmission_rate(house);

	return param("alpha") * member_i.group_rate("community") + hh_roi_inst;
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
		i0 = house.get_infected_status_hhindexes().at(uniform() * house.n_infec_status());
	else
		// amongst all members
		i0 = uniform() * house.size();
	Individual &member_0 = house.get_member(i0);

	// set all to healthy except from introduction
	// assign augmented variables to all
	for (Individual &member : house)
	{
		house.change_status(member, (member == member_0) ? INFECTED_STATUS : HEALTHY_STATUS);
		assign_all_augm(member);
	}

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	// first infection occuring at time t = 0
	const double t0 = member_0.augm("first_pos_date") - member_0.augm("delay_period");
	const double tf = member_0.augm("first_pos_date") + param("study_period");

	// infection of secondary cases following the first case
	for (double t = t0 + 1; (t < tf) && (house.n_healthy_status() > 0); t += 1)
	{
		// store a copy of the vector of infected individuals
		const std::vector<indiv> infectors = house.get_infected_status_hhindexes();

		for (Individual &member_i : house.healthy())
			if (uniform() >= exp(-roi_inst_simulation(infectors, house, member_i, t)))
			{
				member_i.set_augm_value("first_pos_date", t + member_i.augm("delay_period"));
				house.change_status(member_i, INFECTED_STATUS);
			}
	}

	// format infected with unseen onset as healthy
	for (Individual &member : house.infected())
		if (member.augm("first_pos_date") >= tf)
			house.change_status(member, HEALTHY_STATUS);

	// format all augmented variables of susceptibles
	for (Individual &member : house.healthy())
		clear_all_augm(member);
}
