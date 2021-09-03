/**
 * GeoDa TM, Copyright (C) 2011-2015 by Luc Anselin - all rights reserved
 *
 * This file is part of GeoDa.
 *
 * GeoDa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GeoDa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Created: 8/22/2017 lixun910@gmail.com
 */


#ifndef __GEODA_CENTER_JOINCOUNT_RATIO_H__
#define __GEODA_CENTER_JOINCOUNT_RATIO_H__

#include <vector>
#include "../weights/GeodaWeight.h"
#include "../gda_clustering.h"


std::vector<JoinCountRatio> joincount_ratio(const std::vector<int>& clusters, GeoDaWeight* w);

JoinCountRatio all_joincount_ratio(const std::vector<JoinCountRatio>& items);

#endif
