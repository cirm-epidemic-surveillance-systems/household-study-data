/*
 * model.h
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  * *** * E.G.AUGMENTED ONSET - SINGLE INTRODUCTION * *** *  ///////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HOUSEPI__MODEL__H
#define HOUSEPI__MODEL__H

#include "mcmc.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Model : class //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
public:
	//===== Model dependent methods =====//
	double transmission_rate(Household &house) const;
	double roi_inst_simulation(const std::vector<indiv> &infectors,
							   Household &house, indiv i, double t_infec_i) const;

	//===== Override General Methods =====//
	void init_model_for_inference() override final;
	void init_model_for_simulation() override final;
	double compute_house_obs_loglike(Household &house) const override final;
	void simulate_single_house(Household &house) override final;

	//----- Rate Choices -----//
	void assign_group_rate_to_individual(Household &house, indiv i) override final;
	void assign_contact_rate_to_household(Household &house, indiv from, indiv to) override final;
};

#endif // HOUSEPI__MODEL__H
