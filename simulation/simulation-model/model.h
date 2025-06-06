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
							   Household &house,
							   Individual &member_i,
							   double t_infec_i) const;

	//===== Override General Methods =====//
	void init_model_for_inference() override final;
	void init_model_for_simulation() override final;
	void simulate_single_house(Household &house) override final;

	//----- Rate Choices -----//
	void assign_group_rate_to_individual(Individual &member) override final;
	void assign_contact_rate_to_household(Household &house, indiv from, indiv to) override final;
};

#endif // HOUSEPI__MODEL__H
