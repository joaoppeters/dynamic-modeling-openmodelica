package proposta_mit
  model regime_permanente
    //
    //constants & libraries
    //
    import Modelica.Constants.pi;
    import SI = Modelica.SIunits;
    //
    //parameters
    //
    parameter SI.Frequency fb = 60 "base linear frequency";
    parameter SI.AngularFrequency wb = 2 * pi * fb "base angular frequency";
    //
    parameter SI.Power Pn = 50 * 746 "3P induction machine nominal power";
    parameter SI.Voltage Vn = 460 / sqrt(3) "1P induction machine nominal voltage";
    parameter SI.Inertia Jr = 1.662 "rotor mechanical inertia";
    parameter SI.Inertia Js = 15 "blower inertia";
    parameter SI.Inertia J = Jr + Js "total inertia";
    parameter SI.AngularFrequency wrmn = 1720 * 2 * pi / 60 "rotor mechanical nominal speed";
    parameter Integer p = 4 "n# poles";
    //
    parameter SI.Resistance rs = 0.087 "stator winding resistance";
    parameter SI.Resistance rr_ = 0.228 "rotor winding resistance - stator referenced";
    parameter SI.Inductance Lls = 0.302 / wb "stator winding leakage reactance";
    parameter SI.Inductance Llr_ = 0.302 / wb "rotor winding leakage reactance - stator referenced";
    parameter SI.Inductance Lm = 13.08 / wb "magnetizing reactance - stator referenced";
    //
    //variables
    //
    SI.Voltage vqs "stator voltage q-axis";
    SI.Voltage vds "stator voltage d-axis";
    SI.Voltage vqr_ "rotor voltage q-axis";
    SI.Voltage vdr_ "rotor voltage d-axis";
    //
    SI.Current iqs "stator current q-axis";
    SI.Current ids "stator current d-axis";
    SI.Current iqr_ "rotor current q-axis";
    SI.Current idr_ "rotor current d-axis";
    //
    SI.MagneticFlux fqs "stator magnetic flux q-axis";
    SI.MagneticFlux fds "stator magnetic flux d-axis";
    SI.MagneticFlux fqr_ "rotor magnetic flux q-axis";
    SI.MagneticFlux fdr_ "rotor magnetic flux d-axis";
    //
    SI.AngularFrequency w "qd0 referential speed";
    SI.AngularFrequency wr "rotor electrical speed";
    SI.AngularFrequency wrm "rotor mechanical speed";
    //
    SI.Torque Te "electric conjugate";
    SI.Torque Tm "mechanical conjugate";
    //
    SI.Power Pe "electric active power";
    SI.ReactivePower Qe "electric reactive power";
    //
    SI.MagneticFlux fm "magnetic flux";
    SI.MagneticFlux fmq "q-axis magnetic flux";
    SI.MagneticFlux fmd "d-axis magnetic flux";
    //
    SI.Angle theta;
    SI.Angle var_m;
    //
    SI.Voltage va "stator voltage phase A";
    SI.Voltage vb "stator voltage phase B";
    SI.Voltage vc "stator voltage phase C";
    //
    SI.Current ia "stator current phase A";
    SI.Current ib "stator current phase B";
    SI.Current ic "stator current phase C";
    //
    //fim
    //
  initial equation
    der(fqs) = 0;
    der(fds) = 0;
    der(fqr_) = 0;
    der(fdr_) = 0;
    der(wrm) = 0;
  //
  equation
//
//steady-state conditionals
//
    w = wb;
    theta = w * time;
    va = sqrt(2) * Vn * cos(wb * time);
    vb = sqrt(2) * Vn * cos(wb * time - 2 * pi / 3);
    vc = sqrt(2) * Vn * cos(wb * time + 2 * pi / 3);
    ia = iqs * cos(theta) + ids * sin(theta);
    ib = iqs * cos(theta - 2 * pi / 3) + ids * sin(theta - 2 * pi / 3);
    ic = iqs * cos(theta + 2 * pi / 3) + ids * sin(theta + 2 * pi / 3);
    vqs = (2 / 3) * (va * cos(theta) + vb * cos(theta - 2 * pi / 3) + vc * cos(theta + 2 * pi / 3));
    vds = (2 / 3) * (va * sin(theta) + vb * sin(theta - 2 * pi / 3) + vc * sin(theta + 2 * pi / 3));
    vqr_ = 0;
    vdr_ = 0;
//
//mechanical equations
//
    J * der(wrm) = Te - Tm;
    wrm = wrmn;
    wr = p / 2 * wrm;
    Te = 3 * p * (fds * iqs - fqs * ids) / 4;
    Pe = 3 / 2 * (vqs * iqs + vds * ids);
    Qe = 3 / 2 * (vqs * ids - vds * iqs);
//
//electromagnetical equations
//
//stator
    vqs = rs * iqs + w * fds + der(fqs);
    vds = rs * ids - w * fqs + der(fds);
    fqs = Lls * iqs + Lm * (iqs + iqr_);
    fds = Lls * ids + Lm * (ids + idr_);
//rotor
    vqr_ = rr_ * iqr_ + (w - wr) * fdr_ + der(fqr_);
    vdr_ = rr_ * idr_ - (w - wr) * fqr_ + der(fdr_);
    fqr_ = Llr_ * iqr_ + Lm * (iqs + iqr_);
    fdr_ = Llr_ * idr_ + Lm * (ids + idr_);  
  //
    fmq = Lm * (iqs + iqr_);
    fmd = Lm * (ids + idr_);
    fm = sqrt(fmq^2 + fmd^2);
    var_m = -atan(fdr_/fqr_);
    annotation(
      experiment(StartTime = 0, StopTime = 50, Tolerance = 1e-06, Interval = 0.001));
  end regime_permanente;

  model curva_conjugado
    //
    //constants & libraries
    //
    import Modelica.Constants.pi;
    import SI = Modelica.SIunits;
    //
    //parameters
    //
    parameter SI.Frequency fb = 60 "base linear frequency";
    parameter SI.AngularFrequency wb = 2 * pi * fb "base angular frequency";
    //
    parameter SI.Power Pn = 50 * 746 "3P induction machine nominal power";
    parameter SI.Voltage Vn = 460 / sqrt(3) "1P induction machine nominal voltage";
    parameter SI.Inertia Jr = 1.662 "rotor mechanical inertia";
    parameter SI.Inertia Js = 15 "blower inertia";
    parameter SI.Inertia J = Jr + Js "total inertia";
    parameter SI.AngularFrequency wrmn = 1720 * 2 * pi / 60 "rotor mechanical nominal speed";
    parameter Integer p = 4 "n# poles";
    //
    parameter SI.Resistance rs = 0.087 "stator winding resistance";
    parameter SI.Resistance rr_ = 0.228 "rotor winding resistance - stator referenced";
    parameter SI.Inductance Lls = 0.302 / wb "stator winding leakage reactance";
    parameter SI.Inductance Llr_ = 0.302 / wb "rotor winding leakage reactance - stator referenced";
    parameter SI.Inductance Lm = 13.08 / wb "magnetizing reactance - stator referenced";
    //
    //variables
    //
    SI.Voltage vqs "stator voltage quadratic axis";
    SI.Voltage vds "stator voltage direct axis";
    SI.Voltage vqr_ "rotor voltage quadratic axis";
    SI.Voltage vdr_ "rotor voltage direct axis";
    //
    SI.Current iqs "stator current quadratic axis";
    SI.Current ids "stator current direct axis";
    SI.Current iqr_ "rotor current quadratic axis";
    SI.Current idr_ "rotor current direct axis";
    //
    SI.MagneticFlux fqs "stator magnetic flux quadratic axis";
    SI.MagneticFlux fds "stator magnetic flux direct axis";
    SI.MagneticFlux fqr_ "rotor magnetic flux quadratic axis";
    SI.MagneticFlux fdr_ "rotor magnetic flux direct axis";
    //
    SI.AngularFrequency w "qd0 referential speed";
    SI.AngularFrequency wr "rotor electrical speed";
    SI.AngularFrequency wrm "rotor mechanical speed";
    //
    SI.Torque Te "electric conjugate";
    SI.Torque Tm "mechanical conjugate";
    SI.Torque Tn "nominal conjugate";
    parameter SI.Torque Top(start=199.811, fixed=true);
    //
    SI.Power Pe "electric power";
    //
    Modelica.Blocks.Sources.Ramp rampa(duration = 500, height = 1800 * 2 * pi / 60) annotation(
      Placement(visible = true, transformation(origin = {-24, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    //
    //fim
    //
  initial equation
    der(fqs) = 0;
    der(fds) = 0;
    der(fqr_) = 0;
    der(fdr_) = 0;
  equation
//
//steady-state conditionals
//
    Tn = 199.811;
    w = wb;
    vqs = sqrt(2) * Vn;
    vds = 0;
    vqr_ = 0;
    vdr_ = 0;
//
//mechanical equations
//
    J * der(wrm) = Te - Tm;
    wrm = rampa.y;
    wr = p / 2 * wrm;
    Te = 3 * p * (fds * iqs - fqs * ids) / 4;
    Pe = 3 / 2 * (vqs * iqs + vds * ids);
//
//electromagnetical equations
//
//stator
    vqs = rs * iqs + w * fds + der(fqs);
    vds = rs * ids - w * fqs + der(fds);
    fqs = Lls * iqs + Lm * (iqs + iqr_);
    fds = Lls * ids + Lm * (ids + idr_);
//rotor
    vqr_ = rr_ * iqr_ + (w - wr) * fdr_ + der(fqr_);
    vdr_ = rr_ * idr_ - (w - wr) * fqr_ + der(fdr_);
    fqr_ = Llr_ * iqr_ + Lm * (iqs + iqr_);
    fdr_ = Llr_ * idr_ + Lm * (ids + idr_);
    annotation(
      experiment(StartTime = 0, StopTime = 500, Tolerance = 1e-06, Interval = 0.1));
  end curva_conjugado;

  model partida_motor
    //
    //constants & libraries
    //
    import Modelica.Constants.pi;
    import SI = Modelica.SIunits;
    //
    //parameters
    //
    parameter SI.Frequency fb = 60 "base linear frequency";
    parameter SI.AngularFrequency wb = 2 * pi * fb "base angular frequency";
    //
    parameter SI.Power Pn = 50 * 746 "3P induction machine nominal power";
    parameter SI.Voltage Vn = 460 / sqrt(3) "1P induction machine nominal voltage";
    parameter SI.Inertia Jr = 1.662 "rotor mechanical inertia";
    parameter SI.Inertia Js = 15 "blower inertia";
    parameter SI.Inertia J = Jr + Js "total inertia";
    parameter SI.AngularFrequency wrmn = 1720 * 2 * pi / 60 "rotor mechanical nominal speed";
    parameter Integer p = 4 "n# poles";
    //
    parameter SI.Resistance rs = 0.087 "stator winding resistance";
    parameter SI.Resistance rr_ = 0.228 "rotor winding resistance - stator referenced";
    parameter SI.Inductance Lls = 0.302 / wb "stator winding leakage reactance";
    parameter SI.Inductance Llr_ = 0.302 / wb "rotor winding leakage reactance - stator referenced";
    parameter SI.Inductance Lm = 13.08 / wb "magnetizing reactance - stator referenced";
    //
    //variables
    //
    SI.Voltage vqs "stator voltage quadratic axis";
    SI.Voltage vds "stator voltage direct axis";
    SI.Voltage vqr_ "rotor voltage quadratic axis";
    SI.Voltage vdr_ "rotor voltage direct axis";
    //
    SI.Current iqs "stator current quadratic axis";
    SI.Current ids "stator current direct axis";
    SI.Current iqr_ "rotor current quadratic axis";
    SI.Current idr_ "rotor current direct axis";
    //
    SI.MagneticFlux fqs "stator magnetic flux quadratic axis";
    SI.MagneticFlux fds "stator magnetic flux direct axis";
    SI.MagneticFlux fqr_ "rotor magnetic flux quadratic axis";
    SI.MagneticFlux fdr_ "rotor magnetic flux direct axis";
    //
    SI.Angle theta;
    //
    SI.AngularFrequency w "qd0 referential speed";
    SI.AngularFrequency wr "rotor electrical speed";
    SI.AngularFrequency wrm "rotor mechanical speed";
    //
    SI.Torque Te "electric conjugate";
    SI.Torque Tm "mechanical conjugate";
    SI.Torque Tn "nominal conjugate";
    //
    SI.Power Pe "active electric power";
    SI.ReactivePower Qe "reactive electric power";
    //
    SI.Voltage va;
    SI.Voltage vb;
    SI.Voltage vc;
    //
    SI.Current ia;
    SI.Current ib;
    SI.Current ic;
    //
    SI.MagneticFlux fmq;
    SI.MagneticFlux fmd;
    SI.MagneticFlux fm;
    //
    //fim
    //
    //initial equation
    //  der(fqs) = 0;
    //  der(fds) = 0;
    //  der(fqr_) = 0;
    //  der(fdr_) = 0;
    //  der(wrm) = 0;
    //
  equation
  //
  //steady-state conditionals
  //
    w = 0;
    theta = w * time;
    va = sqrt(2) * Vn * cos(wb * time);
    vb = sqrt(2) * Vn * cos(wb * time - 2 * pi / 3);
    vc = sqrt(2) * Vn * cos(wb * time + 2 * pi / 3);
    ia = iqs * cos(theta) + ids * sin(theta);
    ib = iqs * cos(theta - 2 * pi / 3) + ids * sin(theta - 2 * pi / 3);
    ic = iqs * cos(theta + 2 * pi / 3) + ids * sin(theta + 2 * pi / 3);
    vqs = (2 / 3) * (va * cos(theta) + vb * cos(theta - 2 * pi / 3) + vc * cos(theta + 2 * pi / 3));
    vds = (2 / 3) * (va * sin(theta) + vb * sin(theta - 2 * pi / 3) + vc * sin(theta + 2 * pi / 3));
    vqr_ = 0;
    vdr_ = 0;
  //
  //mechanical equations
  //
    Tn = 199.811;
    J * der(wrm) = Te - Tm;
    Tm = Tn*(wrm/wrmn)^2;
    wr = p / 2 * wrm;
    Te = 3 * p * (fds * iqs - fqs * ids) / 4;
    Pe = 3 / 2 * (vqs * iqs + vds * ids);
    Qe = 3 / 2 * (vqs * ids - vds * iqs);
  //
  //electromagnetical equations
  //
  //stator
    vqs = rs * iqs + w * fds + der(fqs);
    vds = rs * ids - w * fqs + der(fds);
    fqs = Lls * iqs + Lm * (iqs + iqr_);
    fds = Lls * ids + Lm * (ids + idr_);
  //rotor
    vqr_ = rr_ * iqr_ + (w - wr) * fdr_ + der(fqr_);
    vdr_ = rr_ * idr_ - (w - wr) * fqr_ + der(fdr_);
    fqr_ = Llr_ * iqr_ + Lm * (iqs + iqr_);
    fdr_ = Llr_ * idr_ + Lm * (ids + idr_);  
  //
    fmq = Lm * (iqs + iqr_);
    fmd = Lm * (ids + idr_);
    fm = sqrt(fmq^2 + fmd^2);
    annotation(
      experiment(StartTime = 0, StopTime = 16, Tolerance = 1e-06, Interval = 0.0005));
  end partida_motor;
  
  model orientacao_campo
    //
    //constants & libraries
    //
    import Modelica.Constants.pi;
    import SI = Modelica.SIunits;
    //
    //parameters
    //
    parameter SI.Frequency fb = 60 "base linear frequency";
    parameter SI.AngularFrequency wb = 2 * pi * fb "base angular frequency";
    //
    parameter SI.Power Pn = 50 * 746 "3P induction machine nominal power";
    parameter SI.Voltage Vn = 460 / sqrt(3) "1P induction machine nominal voltage";
    parameter SI.Inertia Jr = 1.662 "rotor mechanical inertia";
    parameter SI.Inertia Js = 15 "blower inertia";
    parameter SI.Inertia J = Jr + Js "total inertia";
    parameter SI.AngularFrequency wrmn = 1720 * 2 * pi / 60 "rotor mechanical nominal speed";
    parameter Integer p = 4 "n# poles";
    //
    parameter SI.Resistance rs = 0.087 "stator winding resistance";
    parameter SI.Resistance rr_ = 0.228 "rotor winding resistance - stator referenced";
    parameter SI.Inductance Lls = 0.302 / wb "stator winding leakage reactance";
    parameter SI.Inductance Llr_ = 0.302 / wb "rotor winding leakage reactance - stator referenced";
    parameter SI.Inductance Lm = 13.08 / wb "magnetizing reactance - stator referenced";
    //
    //variables
    //
    SI.Voltage vqs "stator voltage quadratic axis";
    SI.Voltage vds "stator voltage direct axis";
    SI.Voltage vqr_ "rotor voltage quadratic axis";
    SI.Voltage vdr_ "rotor voltage direct axis";
    //
    SI.Current iqs "stator current quadratic axis";
    SI.Current ids "stator current direct axis";
    SI.Current iqr_ "rotor current quadratic axis";
    SI.Current idr_ "rotor current direct axis";
    //
    SI.MagneticFlux fqs "stator magnetic flux quadratic axis";
    SI.MagneticFlux fds "stator magnetic flux direct axis";
    SI.MagneticFlux fqr_ "rotor magnetic flux quadratic axis";
    SI.MagneticFlux fdr_ "rotor magnetic flux direct axis";
    //
    SI.Angle theta;
    //
    SI.AngularFrequency w "qd0 referential speed";
    SI.AngularFrequency wr "rotor electrical speed";
    SI.AngularFrequency wrm "rotor mechanical speed";
    //
    SI.Torque Te "electric conjugate";
    SI.Torque Tm "mechanical conjugate";
    SI.Torque Tn "nominal conjugate";
    //
    SI.Power Pe "active electric power";
    SI.ReactivePower Qe "reactive electric power";
    //
    SI.MagneticFlux fmq;
    SI.MagneticFlux fmd;
    SI.MagneticFlux fm;
    //
    SI.Voltage va;
    SI.Voltage vb;
    SI.Voltage vc;
    //
    SI.Current ia;
    SI.Current ib;
    SI.Current ic;
    //
    SI.Current iqs_foc "q-axis foc stator current";
    SI.Current ids_foc "d_axis foc stator current";
    SI.Current iqr_foc "q-axis foc rotor current";
    //
    SI.MagneticFlux fdr_foc "d-axis foc rotor magnetic flux";
    //
    SI.Torque Te_foc "foc electric conjugate";
    //
    //
    //fim
    //
  equation
  //
  //steady-state conditionals
  //
    theta = w * time;
    va = sqrt(2) * Vn * cos(wb * time);
    vb = sqrt(2) * Vn * cos(wb * time - 2 * pi / 3);
    vc = sqrt(2) * Vn * cos(wb * time + 2 * pi / 3);
    ia = iqs * cos(theta) + ids * sin(theta);
    ib = iqs * cos(theta - 2 * pi / 3) + ids * sin(theta - 2 * pi / 3);
    ic = iqs * cos(theta + 2 * pi / 3) + ids * sin(theta + 2 * pi / 3);
    vqr_ = 0;
    vdr_ = 0;
  //
  //mechanical equations
  //
    Tn = 199.811;
    J * der(wrm) = Te - Tm;
    Tm = Tn*(wrm/wrmn)^2;
    wr = p / 2 * wrm;
    Te = 3 * p * (fds * iqs - fqs * ids) / 4;
    Pe = 3 / 2 * (vqs * iqs + vds * ids);
    Qe = 3 / 2 * (vqs * ids - vds * iqs);
  //
  //electromagnetical equations
  //
  //stator
    vqs = rs * iqs + w * fds + der(fqs);
    vds = rs * ids - w * fqs + der(fds);
    fqs = Lls * iqs + Lm * (iqs + iqr_);
    fds = Lls * ids + Lm * (ids + idr_);
  //rotor
    vqr_ = rr_ * iqr_ + (w - wr) * fdr_ + der(fqr_);
    vdr_ = rr_ * idr_ - (w - wr) * fqr_ + der(fdr_);
    fqr_ = Llr_ * iqr_ + Lm * (iqs + iqr_);
    fdr_ = Llr_ * idr_ + Lm * (ids + idr_);  
  //
    fmq = Lm * (iqs + iqr_);
    fmd = Lm * (ids + idr_);
    fm = sqrt(fmq^2 + fmd^2);
  //
  //FOC
  //
    Te_foc = Tn;
    fdr_foc = sqrt((4*Te_foc*rr_)/(3*p*(wb-wrmn*(p/2)))); 
    iqr_foc = -(4*Te_foc)/(fdr_foc*3*p);
    iqs_foc = -iqr_foc*(Llr_+Lm)/Lm;
    ids_foc = fdr_foc/Lm;
    w = wr+(rr_/(Llr_+Lm))*(iqs_foc/ids_foc);
  //
    iqs = iqs_foc;
    ids = ids_foc;
    annotation(
      experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-06, Interval = 0.0001));
  end orientacao_campo;
  annotation(
    uses(Modelica(version = "3.2.3")));
end proposta_mit;