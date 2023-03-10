NEURON {
    SUFFIX dIN_na
    USEION na READ ena WRITE ina
    RANGE gmax, ina, h, m
    RANGE alpha_A_m, alpha_B_m, alpha_C_m, alpha_D_m, alpha_E_m
    RANGE beta_A_m, beta_B_m, beta_C_m, beta_D_m, beta_E_m
    RANGE alpha_A_h, alpha_B_h, alpha_C_h, alpha_D_h, alpha_E_h
    RANGE beta_A_h, beta_B_h, beta_C_h, beta_D_h, beta_E_h
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}

PARAMETER {
    gmax = 0  (S/cm2)

    alpha_A_m = 8.67
    alpha_B_m = 0.0
    alpha_C_m = 1.0
    alpha_D_m = -1.01
    alpha_E_m = -12.56

    beta_A_m = 3.82
    beta_B_m = 0.0
    beta_C_m = 1.0
    beta_D_m = 9.01
    beta_E_m = 9.69

    alpha_A_h = 0.08
    alpha_B_h = 0.0
    alpha_C_h = 0.0
    alpha_D_h = 38.88
    alpha_E_h = 26.0

    beta_A_h = 4.08
    beta_B_h = 0.0
    beta_C_h = 1.0
    beta_D_h = -5.09
    beta_E_h = -10.21
}

ASSIGNED {
    v (mV)
    ena (mV)
    ina (mA/cm2)
    na_minf
    na_hinf
    na_mtau (ms)
    na_htau (ms)
    rate_h
    rate_m
    
}

STATE {
    h 
    m
    
}

INITIAL {
    ena = 50.0
    rates()
    h = na_hinf
    m = na_minf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gmax * pow(m,3) * h * (v - ena)
}

DERIVATIVE states {
    rates()
    m' =  ( na_minf  -  m ) /  na_mtau
    h' =  ( na_hinf - h ) / na_htau
}

UNITSOFF
PROCEDURE rates() {LOCAL alpha_m, beta_m, alpha_h, beta_h
    :alpha_m = alphabeta(8.67,0.0,1.0,-1.01,-12.56)
    :beta_m =  alphabeta(3.82,0.0,1.0,9.01,9.69)
    :alpha_h = alphabeta(0.08,0.0,0.0,38.88,26.0)
    :beta_h =  alphabeta(4.08,0.0,1.0,-5.09,-10.21)

    alpha_m = (alpha_A_m+alpha_B_m*v)/(alpha_C_m+exp((alpha_D_m+v)/alpha_E_m))
    beta_m = (beta_A_m+beta_B_m*v)/(beta_C_m+exp((beta_D_m+v)/beta_E_m))

    alpha_h = (alpha_A_h+alpha_B_h*v)/(alpha_C_h+exp((alpha_D_h+v)/alpha_E_h))
    beta_h = (beta_A_h+beta_B_h*v)/(beta_C_h+exp((beta_D_h+v)/beta_E_h))

    na_mtau = 1/( alpha_m + beta_m)
    na_htau = 1/( alpha_h + beta_h)
    na_minf = alpha_m * na_mtau
    na_hinf = alpha_h * na_htau
}

FUNCTION alphabeta(A,B,C,D,E){
    alphabeta = (A + B*v)/(C + exp((v+D)/E))
}
UNITSON

