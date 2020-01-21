/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
 * GMSHLine.h
 *
 *  Created on: Mar 22, 2012
 *      Author: fischeth
 */

#ifndef GMSHLINE_H_
#define GMSHLINE_H_

#include <ostream>

namespace FileIO {

class GMSHLine {
public:
	GMSHLine(size_t start_point_id, size_t end_point_id);
	virtual ~GMSHLine();
	void write(std::ostream &os, size_t id) const;
	void resetLineData(size_t start_point_id, size_t end_point_id);

private:
	size_t _start_pnt_id;
	size_t _end_pnt_id;
};

}

#endif /* GMSHLINE_H_ */
