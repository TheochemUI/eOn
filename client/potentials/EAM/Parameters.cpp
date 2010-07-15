#include "Parameters.h"

namespace parameters_data {

    element_parameters ep_al = {3.63262, 1.49303,  2.12396,
                                3.45762, 3.54208, 5.555,
                                {-0.00311395, 0.62803, 2.23414, -2.7872, 1.46986,
                                 -0.228432, -0.113191, 0.0556286, -0.0072045}};

    element_parameters ep_au = {0.669134, 1.89707, 2.57021,
                                3.69785, 3.69674, 5.5155,
                                {-0.0073058, -1.02044, 8.23475, -21.8339,
                                 48.3759, -76.0976, 72.5151, -37.1929, 
                                 7.87845}};

    element_parameters ep_pd = {1.6054, 1.55083, 2.35604,
                                3.35148, 3.34715, 5.4120,
                                {0.00364802, 0.177714, -0.0135482,0.634835,
                                 -0.678086, 0.34920, -0.0991826, 0.0148968,
                                 -0.000925276}};
}

element_parameters get_element_parameters(int atomic_number)
{
    using namespace parameters_data;
    switch(atomic_number)
    {
        case 13:
            return ep_al;
        case 46:
            return ep_pd;
        case 79:
            return ep_au;
    }
}
