function rmm11!(dx, x, parameters, t)
    # Reaction: A <-> P
    # Formulation: Reversible Michaelis-Menten
    # Parameters: [1]=Keq [2]=KmA [3]=Vmax [4]=Enzyme/protein amount

    a, p = x
    Keq, KmA, KmP, Vf, enz = parameters

    V = ((a - p / Keq) * enz * Vf / KmA) / (1 + a / KmA + p / KmP)

    dx .= [-V; V]
end

function rhe22!(dx, x, parameters, t)
    # Reaction: A + B <-> P + Q
    # Formulation: Reversible Hill Equation
    # Parameters: [1]=Keq [2]=KmA [3]=KmB [4]=KmP [5]=KmQ [5]=Vmax [6]=Enzyme/protein amount

    a, b, p, q = x
    Keq, KmA, KmB, KmP, KmQ, Vf, enz = parameters

    V = (enz*Vf*(a*b-p*q/Keq)/(KmA*KmB))/ ((1 + a / KmA + p / KmP) * (1 + b / KmB + q / KmQ))

    dx .= [-V; -V; V; V]
end

function rhe22_2enz!(dx, x, Keq1, KmA1, KmB1, KmI1, KmQ1, Vmax1, enz1, Keq2, KmI2, KmP2, KmR2, KmC2, Vmax2, enz2)
    # Reaction 1: A + B <-> I + Q
    # Reaction 2: I + C <-> P + R
    # Formulation: Reversible Hill Equation
    # x and dx order: [A B C I P Q R]
    a, b, c, i, p, q, r = x

    V1 = (enz1*Vmax1*(a*b-i*q/Keq1)/(KmA1*KmB1))/((1+a/KmA1+i/KmI1)*(1+b/KmB1+q/KmQ1))
    V2 = (enz2*Vmax2*(i*c-p*r/Keq2)/(KmI2*KmC2))/((1+i/KmI2+p/KmP2)*(1+c/KmC2+r/KmR2))

    dx .= [-V1; -V1; -V2; V1-V2; V2; V1; V2]
end

function rhe22_2enz_fitEnz1!(dx, x, parameters, t)
    # Reorder the parameters considering that those of enz2 are fixed and those of enz1 are being fit
    #                     FIXED                    |         VARIABLE
    KmI2, KmP2, KmR2, KmC2, Vmax2, enz2, Keq1, Keq2, KmA1, KmB1, KmI1, KmQ1, Vmax1, enz1 = parameters
    rhe22_2enz!(dx, x, Keq1, KmA1, KmB1, KmI1, KmQ1, Vmax1, enz1, Keq2, KmI2, KmP2, KmR2, KmC2, Vmax2, enz2)
end

function rhe22_2enz_fitEnz2!(dx, x, parameters, t)
    # Reorder the parameters considering that those of enz1 are fixed and those of enz2 are being fit
    #                     FIXED                    |         VARIABLE
    KmA1, KmB1, KmI1, KmQ1, Vmax1, enz1, Keq1, Keq2, KmI2, KmP2, KmR2, KmC2, Vmax2, enz2 = parameters
    rhe22_2enz!(dx, x, Keq1, KmA1, KmB1, KmI1, KmQ1, Vmax1, enz1, Keq2, KmI2, KmP2, KmR2, KmC2, Vmax2, enz2)
end

function rhe22_2enz_fitEnz1_2!(dx, x, parameters, t)
    # Reorder the parameters considering that those of enz1 and enz2 are being fit simultaneously
    #      FIXED    |         VARIABLE
    enz2, Keq1, Keq2, KmA1, KmB1, KmI1, KmQ1, Vmax1, enz1, KmI2, KmP2, KmR2, KmC2, Vmax2 = parameters
    rhe22_2enz!(dx, x, Keq1, KmA1, KmB1, KmI1, KmQ1, Vmax1, enz1, Keq2, KmI2, KmP2, KmR2, KmC2, Vmax2, enz2)
end

function cm32!(dx, x, parameters, t)
    # Reaction: A + B + C <-> P + Q
    # Formulation: Common modular
    # Parameters: [1]=Keq [2]=KmA [3]=KmB [4]=KmC [5]=KmP [6]=KmQ [7]=Vmax [8]=Enzyme/protein amount

    # x and dx order: [A B C P Q]
    a, b, c, p, q = x
    Keq, KmA, KmB, KmC, KmP, KmQ, Vf, enz = parameters

    V = (enz*Vf*(a*b*c-p*q/Keq)/(KmA*KmB*KmC))/((1+a/KmA)*(1+b/KmB)*(1+c/KmC) + (1+p/KmP)*(1+q/KmQ) - 1)

    dx .= [-V; -V; -V; V; V]
end

function dm32!(dx, x, parameters, t)
    # Reaction: A + B + C <-> P + Q
    # Formulation: Direct binding modular
    # Parameters: [1]=Keq [2]=KmA [3]=KmB [4]=KmC [5]=KmP [6]=KmQ [7]=Vmax [8]=Enzyme/protein amount

    # x and dx order: [A B C P Q]
    a, b, c, p, q = x
    Keq, KmA, KmB, KmC, KmP, KmQ, Vf, enz = parameters

    V = (enz*Vf*(a*b*c-p*q/Keq)/(KmA*KmB*KmC))/((a*b*c)/(KmA*KmB*KmC) + (p*q)/(KmP*KmQ) + 1)

    dx .= [-V; -V; -V; V; V]
end

function mass_action!(dx, x, parameters, t)
    # Reaction: A -> *
    # Formulation: Consumption of a single metabolite X according to the law of mass action
    # Parameters: [1]=baseline [2]=K [3]=Enzyme/protein amount
    baseline, K = parameters

    dx .= -K*(x .- baseline)
end
