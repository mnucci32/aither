/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (mnucci@pm.me)

    Aither is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Aither is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef LIMITERHEADERDEF
#define LIMITERHEADERDEF

/* This header contains functions to limit the primitive variables during
 * reconstruction
 */

// forward class declarations
class primitive;

primitive LimiterVanAlbada(const primitive &);
primitive LimiterMinmod(const primitive &);
primitive LimiterNone(const int &, const int &);

#endif
