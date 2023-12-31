/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#ifndef GUARD_heterbart_h
#define GUARD_heterbart_h
#include "bart.h"
#include "heterbartfuns.h"
#include "heterbd.h"

class heterbart : public bart
{
public:
   heterbart():bart() {}
   heterbart(size_t m):bart(m) {}
   void pr();
   void draw(double *sigma, rn& gen);
   void draw_conditional(size_t funcP, size_t intvls, double *sigma,
                                    rn& gen, rn& gen_fp, rn& gen_int);
   void draw_conditional_neighborhood(size_t funcP, size_t intvls, double *sigma,
                         rn& gen, rn& gen_fp, rn& gen_int, double sigma_int);
};

#endif
