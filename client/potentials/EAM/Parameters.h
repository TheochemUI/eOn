struct element_parameters {
    double Dm; 
    double alphaM;
    double Rm;
    double beta1;
    double beta2;
    double r_cut;
    double func_coeff[9];
};

element_parameters get_element_parameters(int atomic_number);
