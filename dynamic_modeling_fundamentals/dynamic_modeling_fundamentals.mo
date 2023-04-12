package proposta_fundamentos_conversao
  model pfc_atuador
    import Modelica.Constants.pi;
    import mu0 = Modelica.Constants.mue_0;
    import g = Modelica.Constants.g_n;
    import SI = Modelica.SIunits;
    // sistema eletrico
    SI.Inductance Lh "indutancia equivalente";
    SI.Current i "corrente";
    SI.Force felet "forca eletrica : sistema mecanico";
    SI.Frequency fn=60 "frequencia eletrica";
    SI.AngularFrequency wn=2*pi*fn "frequencia eletrica angular";
    // sistema magnetico
    parameter SI.Area AL = 0.04 * 0.04 "area da perna lateral";
    parameter SI.Area AM = 2 * AL "area da perna central";
    parameter Integer N = 500 "numero de espiras do enrolamento";
    SI.MagneticFlux lambda "fluxo enlacado pela bobina";
    SI.MagneticFlux phi "fluxo magnetico";
    SI.MagneticFluxDensity B "densidade fluxo magnetico";
    // sistema mecanico
    parameter SI.Mass m = 7800 * 0.24 * 0.04 * 0.04 "peso da pecca movel";
    parameter SI.Mass m_l = 0;//581.2592 "peso adicional acoplado pecca movel";
    parameter SI.Mass m_T = m + m_l "massa total do conjunto";
    parameter Real h0 = 4.18879e-3 "distancia do entreferro";
    SI.Velocity v "velocidade de movimento";
    SI.Force fext "forca externa : mecanica"; 
    SI.Force fpeso "forca peso : mecanica";
    SI.Force fpap "forca papelao : entreferro";
    SI.Position y(start = 0);//-5.4214013e-2) "posicao inicial";
    SI.Length h;
  equation
    i = 20;//*sqrt(2)*cos(wn*time);
    lambda = Lh * i;
    felet = lambda ^ 2 / (2* mu0 * AL * N ^ 2);
    fpeso = m_T * g;
    m * der(v) = -fpeso + felet + fext + fpap;
    fext = 0;
    Lh = N ^ 2 * mu0 * AL / h;
    h = h0 - y;
    phi = N * i * mu0 * AL / h;
    B = phi / AM;
  //
    if y>=0 then
      der(y) = 0;
      fpap = -felet + fpeso;
    else
      fpap = 0;
      v = der(y);
    end if;
  //
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
  end pfc_atuador;
end proposta_fundamentos_conversao;