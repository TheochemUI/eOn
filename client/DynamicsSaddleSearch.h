#ifndef DYNAMICSSADDLESEARCH_H
#define DYNAMICSSADDLESEARCH_H

#include "Matter.h"
#include "MinModeSaddleSearch.h"

class DynamicsSaddleSearch
{
    public:
        DynamicsSaddleSearch(Matter *matterPassed, Parameters *parametersPassed);
        ~DynamicsSaddleSearch();

        int run(void);

        Matter *reactant;
        Matter *saddle;
        Matter *product;
        MinModeSaddleSearch *saddleSearch;

    private:
        Parameters *parameters;
};

#endif
