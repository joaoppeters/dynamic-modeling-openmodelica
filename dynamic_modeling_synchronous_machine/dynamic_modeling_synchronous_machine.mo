package proposta_ms
  model maquina_sincrona
    //
    //constants & libraries
    //
    import Modelica.Constants.pi;
    import SI = Modelica.SIunits;
    //
    //parameters
    //
    parameter SI.Frequency fe = 60;
    parameter SI.AngularFrequency we = 2 * pi * fe;
    //
    parameter SI.ApparentPower Sn = 835e6;
    parameter SI.Voltage Vn = 26e3;
    parameter SI.PowerFactor Fpn = 0.85;
    parameter SI.AngularFrequency wrmn = 3600 * 2 * pi / 60;
    parameter Integer p = 2;
    //
    parameter SI.Voltage vfd = 12.2173751750235;
    parameter SI.Inertia J = 0.0658e6;
    parameter Real H = 5.6;
    //
    parameter SI.Resistance rs = 0.00243;
    parameter SI.Inductance Lls = 0.1538/we;
    parameter SI.Resistance rfd = 0.00075;
    parameter SI.Inductance Llfd = 0.1145/we;
    parameter SI.Resistance rkd = 0.01080;
    parameter SI.Inductance Llkd = 0.06577/we;
    parameter SI.Resistance rkq1 = 0.00144;
    parameter SI.Inductance Llkq1 = 0.6578/we;
    parameter SI.Resistance rkq2 = 0.00681;
    parameter SI.Inductance Llkq2 = 0.07602/we;
    parameter SI.Inductance Lmq = Lq - Lls;
    parameter SI.Inductance Lmd = Ld - Lls;
    parameter SI.Inductance Lq = 1.457/we;
    parameter SI.Inductance Ld = 1.457/we;
    //
    //variables
    //
    SI.Voltage vqs;
    SI.Voltage vds;
    SI.Voltage v0s;
    //
    SI.Current iqs;
    SI.Current ids;
    SI.Current i0s;
    SI.Current ifd;
    SI.Current ikd;
    SI.Current ikq1;
    SI.Current ikq2;
    //
    SI.MagneticFlux fqs;
    SI.MagneticFlux fds;
    SI.MagneticFlux f0s;
    SI.MagneticFlux ffd;
    SI.MagneticFlux fkd;
    SI.MagneticFlux fkq1;
    SI.MagneticFlux fkq2;
    //
    SI.Torque Te;
    SI.Torque Tm;
    //
    SI.Voltage vas;
    SI.Voltage vbs;
    SI.Voltage vcs;
    //
    SI.Current ias;
    SI.Current ibs;
    SI.Current ics;
    //
    SI.AngularFrequency wrm;
    SI.AngularFrequency wr;
    //
    SI.Angle delta;
    SI.Angle theta_e;
    SI.Angle theta_r;
    //
    SI.Power Pe;
    SI.ReactivePower Qe;
    //
    //fim
    //
    initial equation
      wr = we;
      der(fqs) = 0;
      der(fds) = 0;
      der(f0s) = 0;
      der(fkq1) = 0;
      der(fkq2) = 0;
      der(ffd) = 0;
      der(fkd) = 0;
      der(wrm) = 0;
    //
    equation
    vas = Vn*sqrt(2)/sqrt(3)*cos(we*time);
    vbs = Vn*sqrt(2)/sqrt(3)*cos(we*time - 2*pi/3);
    vcs = Vn*sqrt(2)/sqrt(3)*cos(we*time + 2*pi/3);
    ias = i0s + iqs*cos(theta_r) + ids*sin(theta_r);
    ibs = i0s + iqs*cos(theta_r-2*pi/3) + ids*sin(theta_r-2*pi/3);
    ics = i0s + iqs*cos(theta_r+2*pi/3) + ids*sin(theta_r+2*pi/3);
    //
    vqs = (2/3)*(vas*cos(theta_r) + vbs*cos(theta_r-2*pi/3) + vcs*cos(theta_r+2*pi/3));
    vds = (2/3)*(vas*sin(theta_r) + vbs*sin(theta_r-2*pi/3) + vcs*sin(theta_r+2*pi/3));
    v0s = (2/3)*(vas + vbs + vcs);
    //
    vqs = -iqs*rs + fds*we + der(fqs);
    vds = -ids*rs - fqs*we + der(fds);
    v0s = -i0s*rs + der(f0s);
    0 = ikq1*rkq1 + der(fkq1);
    0 = ikq2*rkq2 + der(fkq2);
    vfd = ifd*rfd + der(ffd);
    0 = ikd*rkd + der(fkd);
    //
    fqs = -(Lls+Lmq)*iqs + Lmq*(ikq1+ikq2);
    fds = -(Lls+Lmd)*ids + Lmd*(ifd+ikd);
    f0s = Lls*i0s;
    fkq1 = (Llkq1+Lmq)*ikq1 + Lmq*(ikq2-iqs);
    fkq2 = (Llkq2+Lmq)*ikq2 + Lmq*(ikq1-iqs);
    ffd = (Llfd+Lmd)*ifd + Lmd*(-ids+ikd);
    fkd = (Llkd+Lmd)*ikd + Lmd*(-ids+ifd);
    //
    Tm = if time < 1 then 0 else 1.11e6;
    //
    delta = theta_r - theta_e;
    der(theta_e) = we;
    der(theta_r) = wr;
    wr = (p/2)*wrm;
    Te = (3/2)*(p/2)*(iqs*fds - ids*fqs);
    J*der(wrm) = Tm - Te;
    //
    Pe = (3/2)*(ids*vds + iqs*vqs);
    Qe = (3/2)*(ids*vqs - iqs*vds);
    annotation(
      experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.0001));
  end maquina_sincrona;
  annotation(
    uses(Modelica(version = "3.2.3")));
end proposta_ms;