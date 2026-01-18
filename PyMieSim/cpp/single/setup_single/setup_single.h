#pragma once

#include <functional>
#include <memory>
#include <cstdarg>
#include <cstdio>
#include "single/scatterer/sphere/sphere.h"
#include "single/scatterer/cylinder/cylinder.h"
#include "single/scatterer/coreshell/coreshell.h"
#include "single/source/source.h"
#include "single/detector/detector.h"

class Setup
{
    public:
        std::shared_ptr<BaseSource> source;
        std::shared_ptr<BaseScatterer> scatterer;
        std::shared_ptr<Photodiode> detector;
        bool debug_mode = false;
        Setup() = default;
        Setup(
            std::shared_ptr<BaseScatterer> _scatterer,
            std::shared_ptr<BaseSource> _source,
            std::shared_ptr<Photodiode> _detector,
            bool _debug_mode = false)
            : source(_source), scatterer(_scatterer), detector(_detector), debug_mode(_debug_mode)
        {
            scatterer->source = source;
            detector->source = source;
        }

        Setup(
            std::shared_ptr<BaseScatterer> _scatterer,
            std::shared_ptr<BaseSource> _source,
            bool _debug_mode = false)
            : source(_source), scatterer(_scatterer), debug_mode(_debug_mode)
        {
            scatterer->source = source;
            detector->source = source;
        }


        double get_data(std::string data_name) const {
            if (data_name == "coupling") {
                return detector->get_coupling(*scatterer);
            } else if (data_name == "Qsca") {
                return scatterer->get_Qsca();
            } else if (data_name == "Qext") {
                return scatterer->get_Qext();
            } else if (data_name == "Qback") {
                return scatterer->get_Qback();
            } else if (data_name == "Qforward") {
                return scatterer->get_Qforward();
            } else if (data_name == "g") {
                return scatterer->get_g();
            } else {
                throw std::invalid_argument("Unknown data name: " + data_name);
            }
        }



};
