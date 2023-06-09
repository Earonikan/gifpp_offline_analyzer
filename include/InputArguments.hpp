#pragma once

#include "stdafx.hpp"
#include "ParameterParser.hpp"
#include "DataStructs.hpp"

class InputArguments
{
private:
    Parameters parameters;
public:
    InputArguments(int argv, char *argc[]);
    Parameters Get() {return parameters;}
};