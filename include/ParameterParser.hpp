#pragma once

#include "stdafx.hpp"


class ParameterParser
{
private:
    std::map<std::string, int> data_map;

public:
    ParameterParser(std::vector<std::string> vec_str)
    {
        int k = 1;
        for (auto i : vec_str) data_map[i] = k++;
    }
    int ParseThis(std::string input)
    {
        return data_map[input];
    }


};