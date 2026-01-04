
# Sami2 is Another Model of the Ionosphere (SAMI2): A new low-latitude ionosphere model

# J. D. Huba and G. Joyce

Plasma Physics Division, Naval Research Laboratory, Washington, D.C.

#### J. A. Fedder

Leading Edge Technology, Inc., Washington, D.C.

**Abstract.** A new low-latitude ionospheric model has been developed at the Naval Research Laboratory: Sami2 is Another Model of the Ionosphere (SAMI2). SAMI2 treats the dynamic plasma and chemical evolution of seven ion species  $(H^+, He^+, N^+, O^+, N_2^+, NO^+, and O_2^+)$  in the altitude range  $\sim 100$  km to several thousand kilometers. The ion continuity and momentum equations are solved for all seven species; the temperature equation is solved for H<sup>+</sup>, He<sup>+</sup>, O<sup>+</sup>, and the electrons. SAMI2 models the plasma along the Earth's dipole field from hemisphere to hemisphere, includes the  $\mathbf{E} \times \mathbf{B}$  drift of a flux tube (both in altitude and in longitude), and includes ion inertia in the ion momentum equation for motion along the dipole field line. The final point is relevant for plasma dynamics at very high altitudes where ion inertia can be important. For example, we have found that ion sound waves, which are supported by ion inertia, may be generated in the topside ionosphere (> 1000 km) at sunrise and sunset [Huba et al., 2000b]. The neutral species are specified using the Mass Spectrometer Incoherent Scatter model (MSIS86) and the Horizontal Wind Model (HWM93). In this paper we describe in detail the SAMI2 model and present representative results from the model.

### 1. Introduction

Over the past two decades a number of computational models of the ionosphere have been developed. An excellent overview of the most widely used models is given in STEP: Handbook of Ionospheric Models [Schunk, 1996] and by Anderson et al. [1998]. In general, ionospheric models treat the global ionosphere in three parts: low latitude, midlatitude, and high latitude. Low-latitude models (e.g., the Phillips Laboratory global theoretical ionosphere model (GTIM) [Anderson, 1971; Anderson et al., 1996, the University of Alabama field line interhemispheric plasma model (FLIP) [Richards and Torr, 1996], and the Sheffield University plasmasphere-ionosphere model (SUPIM) [Bailey and Balan, 1996]) consider the plasma dynamics along an entire field line from hemisphere to hemisphere. Mid-latitude models (e.g., GTIM [Decker et al., 1994] and the Utah State University time-dependent ionosphere model (TDIM) [Schunk, 1988; Schunk and Sojka, 1996]) typically have the upper boundary < 1000 km; additional boundary conditions at the upper boundary (e.g., particle flux, heat flux) must be

imposed that are generally not consistent with interhemispheric transport. Finally, high-latitude models (e.g., TDIM, GTIM) also have an upper boundary typically set at 1000 km. However, an important aspect of high latitude models is that magnetospheric effects need to be included: for example, the magnetospheric electric field and auroral precipitation effects. A common feature of these models is that they use empirical neutral atmosphere models such as the Mass Spectrometer Incoherent Scatter model (MSIS) [Hedin, 1987] and the Horizontal Wind Model (HWM) [Hedin et al., 1991] or observed data to specify neutral atmosphere densities and winds. There are also two global ionospheric models that solve the neutral atmosphere equations to determine the neutral atmosphere composition and dynamics, and self-consistently couple this solution to ionospheric dynamics: the National Center for Atmospheric Research (NCAR) thermosphereionosphere general circulation model (TIGCM) [Roble et al., 1988; Roble, 1996] and coupled thermosphereionosphere model (CTIM) [Fuller-Rowell et al., 1996].

Although several of the low-latitude models (GTIM, FLIP, and SUPIM) use the same basic geometry, an interhemispheric flux tube, there are considerable differences between the models. GTIM solves the ion continuity and momentum equations to determine the O<sup>+</sup> ion density with the ion and electron temperatures specified by an empirical model. For flux tubes with an apex altitude greater than  $\sim 1000$  km, GTIM does not solve along the entire interhemispheric flux tube but rather specifies an O<sup>+</sup> density or flux at the upper boundary. GTIM includes an electric field model which allows the plasma to  $\mathbf{E} \times \mathbf{B}$  drift in the "radial" direction. FLIP solves the ion continuity and momentum equations for H<sup>+</sup>, O<sup>+</sup>, He<sup>+</sup>, and N<sup>+</sup>, as well as the electron and ion thermal equations. H<sup>+</sup> and O<sup>+</sup> are treated as major ions, while He<sup>+</sup> and N<sup>+</sup> are treated as minor ions. FLIP also solves ion continuity and momentum equations to obtain the densities of the minor neutral species  $N(^{2}D)$ ,  $N(^4S)$ , NO, and thermal balance equations for the first six vibrational levels of N<sub>2</sub>. FLIP is able to model interhemispheric flux tubes that extend into the plasmasphere; that is, the apex altitude of a flux tube can be > 10,000 km. Recently, FLIP was upgraded to include the  $\mathbf{E} \times \mathbf{B}$  drift motion of a flux tube. SUPIM solves the continuity and momentum equations for H<sup>+</sup>, O<sup>+</sup>, He<sup>+</sup>,  $N_2^+$ ,  $O_2^+$ , and  $NO^+$  and the ion temperature equations for H<sup>+</sup>, O<sup>+</sup>, and He<sup>+</sup> as well as the electron temperature equation. SUPIM models flux tubes that extend into the plasmasphere, and it also includes the  $\mathbf{E} \times \mathbf{B}$ drift motion of a flux tube.

Recently, a new low-latitude model of the ionosphere has been developed at the Naval Research Laboratory (NRL): Sami2 is Another Model of the Ionosphere (SAMI2). It is an outgrowth of a midlatitude ionospheric model developed at NRL in the mid-1970s [Oran et al., 1974]. SAMI2 treats the dynamic plasma and chemical evolution of seven ion species (H<sup>+</sup>, He<sup>+</sup>,  $N^+$ ,  $O^+$ ,  $N_2^+$ ,  $NO^+$ , and  $O_2^+$ ). The ion continuity and momentum equations are solved for all seven ion species. All ion species are treated the same; there are no major or minor ion species. Thermal balance equations are solved for three ion species (H<sup>+</sup>, He<sup>+</sup>, and O<sup>+</sup>) and for the electrons. The neutral atmosphere is specified using the empirical codes MSIS86 and HWM93. SAMI2 models the plasma along a magnetic flux tube from hemisphere to hemisphere. It includes the  $\mathbf{E} \times \mathbf{B}$ drift of a flux tube (both in altitude and in longitude), and it also includes ion inertia in the ion momentum equations for plasma motion along the dipole field line. The previously discussed ionospheric models assume that the ions are collisional at all altitudes and neglect ion inertia. To our knowledge, SAMI2 is the first lowlatitude ionospheric model to include the ion inertial terms in the momentum equations.

The purpose of this paper is to describe SAMI2 in detail and to provide a set of representative results. SAMI2 solves the continuity, momentum, and temperature equations for ions and electrons in an offset, tilted magnetic dipole coordinate system. The details of this coordinate system are described in section 2. This coordinate system is used for two reasons. First, it is a good approximation to the Earth's geomagnetic field [Bailey and Balan, 1996] and second, the parallel and perpendicular dynamics of the system decouple naturally in this reference frame.

The ion and electron equations solved in SAMI2 are presented in section 3. A time-splitting scheme is used: the dynamics of the system is first solved for motion along the geomagnetic field and then for motion transverse to the geomagnetic field. The motion of the plasma along the geomagnetic field is described by a set of advection/diffusion equations. The advection terms are solved using an implicit donor cell method, while the diffusion terms are backward biased for stability. The motion of the plasma orthogonal to the geomagnetic field is determined from the conservation laws for mass and magnetic flux rather than by directly solving the equations. This method is very efficient but assumes that the plasma only  $\mathbf{E} \times \mathbf{B}$  drifts transverse to the magnetic field. Further details of the numerical scheme are discussed in section 4.

The photodeposition model used in SAMI2 is described in section 5. The daytime photoionization model uses the solar EUV flux model for aeronomic calculations (EUVAC) developed by Richards et al. [1994]. The photoionization and photoabsorption rates are taken from Torr and Torr [1982]. The nighttime photoionization model uses the EUV flux prescribed by Strobel et al. [1974], and the photoionization cross sections are obtained from Oran et al. [1974]. The chemistry model is described in section 6. We include 21 chemical reactions and recombination. The rates are taken from Schunk and Sojka [1996] and Bailey and Balan [1996]. The collisional coefficients (e.g., collision frequencies, thermal conductivities) used in SAMI2 are given in Appendix C. Finally, in section 7 we present simulation results from SAMI2.

# 2. Dipole Coordinate System

Ionospheric data are based on a spherical coordinate system aligned with the rotation axis of the Earth. The coordinates are the geographic longitude  $\phi_q$  (with  $\phi_q = 0^0$  at the meridian passing through Greenwich, England; it is positive moving eastward), the geographic latitude  $\theta_g$  (with  $\theta_g = 0^0$  at the equator,  $\theta_g = 90^0$  at the North Pole, and  $\theta_g = -90^{\circ}$  at the South Pole), and the altitude  $h_g$  (where  $h_g = r_g - R_E$  and  $R_E$  is the radius of the Earth). The subscript g denotes geographic coordinates. However, it is much easier to solve the plasma dynamic equations in the frame of the geomagnetic field, because the plasma motions parallel and perpendicular to the magnetic field are weakly coupled. An appropriate coordinate system is an eccentric dipole system which approximately models the Earth's geomagnetic field. The eccentric dipole system is also referred to as an offset, tilted dipole system. The coordinates of the system are typically denoted by q (the coordinate along the field line), p (the coordinate orthogonal to qin the "radial" direction), and  $\phi_d$  (the magnetic longitude, where the subscript d denotes dipole). The transformation from the spherical geographic coordinate system to the offset, tilted dipole coordinate system is accomplished in three stages: (1) spherical geographic to spherical tilted, (2) spherical tilted to spherical eccentric, and (3) spherical eccentric to dipole. We discuss in detail each of these stages. The first two transformations are described in Appendix A while the last transformation is discussed here.

The dipole coordinate system  $(q, p, \phi)$  is defined in terms of the spherical eccentric coordinates by the fol-

lowing:

$$q = \frac{R_E^2}{r_c^2} \cos \theta_e \ , \tag{1}$$

$$p = \frac{r_e}{R_E} \frac{1}{\sin^2 \theta_e} \,, \tag{2}$$

$$\phi = \phi_e \ , \tag{3}$$

where  $R_E$  is the Earth's radius. This is shown in Figure 1

The reverse transformation, from q and p to  $\theta_e$  and  $r_e$ , requires solving the following equation:

$$q^2 \left(\frac{r_e}{R_E}\right)^4 + \frac{1}{p} \left(\frac{r_e}{R_E}\right) - 1 = 0. \tag{4}$$

This is done numerically using Muller's method.

The magnetic field in this coordinate system is

$$\mathbf{B} = \mathbf{B}_{\mathbf{q}} \, \mathbf{e}_{\mathbf{q}} \,\,, \tag{5}$$

where

$$B_q = B_0 \left( 1 + 3\cos^2 \theta_e \right)^{1/2} \left( \frac{R_E}{r_e} \right)^3 .$$
 (6)

This is the key feature of the dipole coordinate system: the magnetic field is aligned with one of the coordinates.

The grid points on a field line are selected so that the smallest spacing between points is at the lowest altitudes, with the spacing increasing approximately exponentially with altitude. We follow the procedure of Millward et al. [1996], which is as follows. First, we determine the q values at the north  $(q_N)$  and south  $(q_S)$ ends of the field line. The ends of the field line correspond to a predetermined minimum altitude  $r_{\min}$ , usually set to 90 km. The q values are determined from (1) subject to the constraint that p = constant. The value q = 0 defines the magnetic equator with negative values of q in the Southern Hemisphere and positive values of q in the Northern Hemisphere. The equation for q must be solved iteratively to obtain  $q_N$  and  $q_S$ , because  $r_{\min}$ is defined in the geographic coordinate system while qis defined in the eccentric coordinate system.

Following *Millward et al.* [1996] we define a constant stride along the field line from the nonlinear equation

$$x_i = \frac{\sinh \gamma q_i}{\sinh \gamma q_S} \,, \tag{7}$$

where  $x_i = x_{i-1} + c$  and c is a constant. If we want  $N_z$  points on the field line, then

$$c = \frac{1 - \sinh \gamma q_N}{\sinh \gamma q_S} \frac{1}{N_z - 1} \ . \tag{8}$$

Equation (7) must be solved for  $q_i$  for each value of  $x_i$ ; we use a Muller scheme for this solution. The parameter  $\gamma$  is a constant that determines the distribution of points along the field line. Larger values of  $\gamma$  produce a larger separation of points at higher altitudes. In practice, a value of  $\gamma$  between 5 and 7 is satisfactory. Finally, in SAMI2 we use a dimensional variable s along the magnetic field instead of q; we define the s coordinates as  $s_i = q_i R_E$ .

# 3. Basic Equations

In deriving the plasma transport equations for the dipole coordinate system, the curvilinear factors are needed. These are given in Appendix B.

#### 3.1. Ion Continuity Equation

The ion continuity equation for each ion species i is

$$\frac{\partial n_i}{\partial t} + \nabla \cdot (n_i \mathbf{V_i}) = \mathcal{P}_i - \mathcal{L}_i n_i , \qquad (9)$$

which can be written in dipole coordinates as

$$\frac{\partial n_i}{\partial t} + b_s^2 \frac{\partial}{\partial s} \frac{n_i V_{is}}{b_s} + \frac{p^2}{R_E^2} b_s^2 \frac{\partial}{\partial p} \frac{r \sin \theta \, n_i V_{ip}}{B_s} + \frac{1}{r \sin \theta} \frac{\partial n_i V_{i\phi}}{\partial \phi} = \mathcal{P}_i - \mathcal{L}_i n_i \,, \tag{10}$$

where  $\mathcal{P}_i$  are the ion production terms and  $\mathcal{L}_i$  are the ion loss terms. These terms involve photoionization, radiative recombination, and chemistry; they are discussed in detail in sections 4 and 5. The second term on the left-hand side of (10) is the field-aligned transport term and the remaining two terms describe motion perpendicular to the magnetic field. The spherical coordinates r and  $\theta$  refer to the eccentric spherical coordinate system; the subscript e has been dropped for simplicity.

#### 3.2. Ion Momentum Equation

The ion momentum equation is

$$\frac{\partial \mathbf{V_i}}{\partial t} + \mathbf{V_i} \cdot \nabla \mathbf{V_i} = -\frac{1}{\rho_i} \nabla \mathbf{P_i} + \frac{\mathbf{e}}{\mathbf{m_i}} \mathbf{E} + \frac{\mathbf{e}}{\mathbf{m_i} \mathbf{c}} \mathbf{V_i} \times \mathbf{B} + \mathbf{g}$$
$$-\nu_{in} (\mathbf{V_i} - \mathbf{V_n}) - \sum_{\mathbf{j}} \nu_{i\mathbf{j}} (\mathbf{V_i} - \mathbf{V_j}) , \qquad (11)$$

where  $\nu_{in}$  is the ion-neutral collision frequency,  $\nu_{ij}$  is the ion-ion collision frequency, and the summation is over ion species  $j \neq i$ .

The field-aligned (s) and cross-field  $(p \text{ and } \phi)$  components of the ion momentum equation are the following:

$$\frac{\partial V_{is}}{\partial t} + (\mathbf{V}_i \cdot \nabla)V_s = -\frac{1}{n_i m_i} b_s \frac{\partial P_i}{\partial s} + \frac{e}{m_i} E_s 
+ g_s - \nu_{in} (V_{is} - V_{ns}) - \sum_j \nu_{ij} (V_{is} - V_{js}) , \qquad (12)$$

$$\frac{\partial V_{ip}}{\partial t} + (\mathbf{V}_i \cdot \nabla)V_{ip} = -\frac{1}{n_i m_i} \frac{p^2}{R_E^2} b_s r \sin \theta \frac{\partial P_i}{\partial p} + \frac{e}{m_i} E_p 
+ g_p - \frac{eB_s}{m_i c} V_{i\phi} - \nu_{in} (V_{ip} - V_{np}) - \sum_j \nu_{ij} (V_{ip} - V_{jp}) , \qquad (13)$$

$$\frac{\partial V_{i\phi}}{\partial t} + (\mathbf{V}_i \cdot \nabla)V_{\phi} = -\frac{1}{n_i m_i} \frac{1}{r \sin \theta} \frac{\partial P_i}{\partial \phi} + \frac{e}{m_i} E_{\phi} 
+ g_{\phi} + \frac{eB_s}{m_i c} V_{ip} - \nu_{in} (V_{i\phi} - V_{n\phi}) - \sum_j \nu_{ij} (V_{i\phi} - V_{j\phi}) , \qquad (14)$$
where  $\mathbf{B} = B_s \mathbf{e}_s$ ,  $b_s = B_s / B_0$ , and  $\Omega_i = eB_s / m_i c$ .

#### 3.3. Electron Momentum Equation

The electron momentum equation used is

$$0 = -\frac{1}{n_e m_e} \nabla \mathbf{P}_e - \frac{e}{m_e} \mathbf{E} - \frac{e}{m_e c} \mathbf{V}_e \times \mathbf{B} . \tag{15}$$

Electron inertia is neglected because of the small electron mass; electron collisional terms are neglected because  $\nu_e \ll \Omega_e$ , where  $\nu_e$  denotes the electron collision frequencies and  $\Omega_e$  is the electron cyclotron frequency. The components of (15) are

$$0 = -\frac{1}{n_a m_e} b_s \frac{\partial P_e}{\partial s} - \frac{e}{m_e} E_s , \qquad (16)$$

$$0 = -\frac{1}{n_e m_e} \frac{p^2}{R_E^2} b_s r \sin \theta \frac{\partial P_e}{\partial p} - \frac{e}{m_e} E_p + \frac{eB_s}{m_e c} V_{e\phi} , \quad (17)$$

$$0 = -\frac{1}{n_e m_e} \frac{1}{r \sin \theta} \frac{\partial P_e}{\partial \phi} - \frac{e}{m_e} E_{\phi} - \frac{e B_q}{m_e c} V_{ep} . \tag{18}$$

Solving (16) – (18) leads to

$$E_s = -\frac{1}{n_e} b_s \frac{\partial P_e}{\partial s} , \qquad (19)$$

$$V_{ep} = -\frac{cE_{\phi}}{B_s} - \frac{c}{en_e B_s} \frac{1}{r \sin \theta} \frac{\partial P_e}{\partial \phi} , \qquad (20)$$

$$V_{e\phi} = \frac{cE_p}{B_s} + \frac{c}{en_e B_s} \frac{p^2}{R_E^2} b_s r \sin \theta \frac{\partial P_e}{\partial p} \ . \tag{21}$$

The components of  $V_e$ , to lowest order, are simply

$$V_p = \frac{cE_\phi}{B_r} \ , \tag{22}$$

$$V_{\phi} = -\frac{cE_p}{B_s} \,, \tag{23}$$

because the pressure terms in (20) and (21) scale as  $\rho_i/R_E << 1$ , where  $\rho_i = C_s/\Omega_i$  and  $C_s = (T_e/m_i)^{1/2}$ .

#### 3.4. Ion Temperature Equation

We write the temperature equation as [Banks and Kockarts, 1973; Millward et al., 1996]

$$\frac{\partial T_i}{\partial t} + \mathbf{V}_i \cdot \nabla T_i + \frac{2}{3} T_i \nabla \cdot \mathbf{V}_i + \frac{2}{3} \frac{1}{n_i k} \nabla \cdot \mathbf{Q}_i = Q_{in} + Q_{ii} + Q_{ie} ,$$
(24)

where the heating terms are due to ion-neutral collisions  $(Q_{in})$ , ion-ion collisions  $(Q_{ij})$ , and ion-electron collisions  $(Q_{ie})$ .

The heat flux  $\mathbf{Q}_i$  is defined as

$$\mathbf{Q}_i = -\kappa_i \nabla T_i \ , \tag{25}$$

and  $\kappa_i$  is the thermal conductivity (given in Appendix C). The thermal conductivity occurs primarily only along the magnetic field line, i.e.,

$$\nabla \cdot \mathbf{Q}_i = -\nabla \cdot \kappa_i \nabla T_i \simeq -b_s^2 \frac{\partial}{\partial s} \kappa_i \frac{\partial T_i}{\partial s} \ . \tag{26}$$

The ion temperature equation can be written as

$$\frac{\partial T_i}{\partial t} + \frac{p^2}{R_E^2} b_s r \sin \theta \frac{\partial T_i}{\partial p} V_{ip} + b_s V_{is} \frac{\partial T_i}{\partial s} + \frac{1}{r \sin \theta} \frac{\partial T_i}{\partial \phi} V_{i\phi} + \frac{2}{3} T_i \left( \frac{p^2}{R_E^2} b_s^2 \frac{\partial}{\partial p} \frac{r \sin \theta V_{ip}}{b_s} + b_s^2 \frac{\partial}{\partial s} \frac{V_{is}}{b_s} + \frac{1}{r \sin \theta} \frac{\partial V_{i\phi}}{\partial \phi} \right) -$$

$$\frac{2}{3}\frac{1}{n_{ik}}b_{s}^{2}\frac{\partial}{\partial s}\kappa_{i}\frac{\partial T_{i}}{\partial s} = Q_{in} + Q_{ij} + Q_{ie} . \tag{27}$$

The various heating terms used in SAMI2 are given as follows:

$$Q_{in} = \sum_{q} \frac{2m_i m_q}{\left(m_i + m_q\right)^2} \nu_{iq} \times$$

$$\left[k\left(T_{q}-T_{i}\right)+\frac{1}{3}m_{q}\left|\mathbf{V_{q}}-\mathbf{V_{i}}\right|^{2}\right],\qquad(28)$$

$$Q_{ii} = \sum_{i} \frac{2.2 \times 10^{-4} n_j}{A_i A_j \left( T_i / A_i + T_j / A_j \right)^{3/2}} \left( T_j - T_i \right) , \quad (29)$$

$$Q_{ie} = \frac{7.7 \times 10^{-6} n_e}{A_i T_e^{3/2}} \left( T_e - T_i \right) , \qquad (30)$$

where q denotes a summation over neutrals and j denotes a summation over ions. In SAMI2 we solve the temperature equation for three ion species:  $\mathrm{H^+}$ ,  $\mathrm{He^+}$ , and  $\mathrm{O^+}$ . We set the temperature of the molecular ions  $\mathrm{N_2^+}$ ,  $\mathrm{NO^+}$ , and  $\mathrm{O_2^+}$  equal to the  $\mathrm{O^+}$  temperature.

#### 3.5. Electron Temperature Equation

The electron temperature equation is

$$\frac{\partial T_e}{\partial t} - \frac{2}{3} \frac{1}{n_e k} b_s^2 \frac{\partial}{\partial s} \kappa_e \frac{\partial T_e}{\partial s} = Q_{en} + Q_{ei} + Q_{\text{phe}} , \quad (31)$$

where  $Q_{en}$  is heating due to electron-neutral collisions,  $Q_{ei}$  is due to electron-ion collisions,  $Q_{phe}$  is due to photoelectron heating, and we have taken  $V_{es} = 0$ . The expression for the parallel electron thermal conductivity  $\kappa_e$  is given in Appendix C.

The collisional heating terms are [Banks and Kockarts, 1973; Bailey and Balan, 1996; Millward et al., 1996]

$$Q_{en} = \sum_{q} \frac{2m_e m_q}{(m_e + m_q)^2} \nu_{eq} (T_q - T_e) , \qquad (32)$$

$$Q_{ei} = \sum_{j} \frac{7.7 \times 10^{-6} n_j}{A_j T_e^{3/2}} (T_j - T_e) , \qquad (33)$$

where q denotes a summation over neutrals and j denotes a summation over ions.

The photoelectron heating model is taken from *Mill-ward et al.* [1996]. At altitudes less than 300 km it is assumed that the photoelectron energy is deposited locally with a rate given by

$$Q_{\rm phe} = \frac{2}{3n_e} \epsilon P_{\rm photo} , \qquad (34)$$

where  $P_{photo}$  is the total photoionization rate. The efficiency factor  $\epsilon$  is an empirical relationship developed by Swartz and Nisbet [1972]:

$$\epsilon = \exp\left[-p(x)\right] , \qquad (35)$$

where

$$p(x) = 12.75 + 6.941x + 1.66x^{2} + 0.08034x^{3} + 0.001996x^{4},$$
 (36)

$$x = \ln\left(\frac{n_e}{n(O_2) + n(N_2) + 0.1n(O)}\right) . \tag{37}$$

At altitudes greater than 300 km it is assumed that the photoelectrons can propagate along the magnetic field, which leads to nonlocal heating. The expression used for  $Q_{\rm phe}$  is

$$Q_{\rm phe} = \frac{2}{3} \frac{B}{B_{300}} q_{300} \exp\left(-C \int n_e dl\right) , \qquad (38)$$

where  $q_{300}$  is the heating rate per electron at 300 km,  $B_{300}$  is the magnetic field strength at 300 km, and dl is the differential length along the geomagnetic field line. The constant C is taken to be  $3 \times 10^{-14}$  cm<sup>2</sup>. This is a rather simple model that will eventually be replaced by a more sophisticated model [Richards and Torr, 1996].

#### 4. Numerical Methods

Existing low-latitude models of the ionosphere assume that the ions are collisional at all altitudes; the equations are solved in a fully implicit manner along the magnetic field. This is done by neglecting ion inertia in (12) and explicitly solving for the ion velocity  $V_{is}$ . The ion velocity  $V_{is}$  is substituted into the ion continuity equation (10), and an implicit equation for  $n_i$  can be obtained. This method requires expansion of the equations (e.g., usually in a Taylor series), which become long and complex. In practice, the convergence to a stable solution is not always guaranteed and, in fact, can fail when ion temperature advection terms are included. Moreover, it is difficult to modify fully implicit codes to include additional ion species (for instance, to model ionospheric disturbances created by shuttle engine burns or chemical releases).

In contrast, the advantages of the numerical methods used in SAMI2 are simplicity, robustness, and flexibility. The equations are not expanded but solved in a straightforward manner. As long as the time step is sufficiently small, i.e., satisfies the Courant condition  $\Delta t < (\Delta l_i/V_{is})_{min}$ , the code solves the entire set of equations with no difficulty. In addition, it is relatively easy to include additional neutral and ion species. The main disadvantage of SAMI2 is that a relatively small time step must be used. Fully implicit codes can use time steps in the range 5-30 min. For altitudes less than  $\sim 4000$  km, SAMI2 uses  $\sim 201$  grid points for adequate resolution and requires a time step of  $\sim 1-12$  s. However, flux tubes that reach higher altitudes require more grid points and necessitate an even smaller time step.

SAMI2 solves the system of equations using a time-splitting method. Schematically, the equations are solved in two steps: (1) parallel transport:  $t_0 + \Delta t \rightarrow t^*$  and (2) perpendicular transport:  $t^* + \Delta t \rightarrow t_1$ . The density, velocity, and temperature are advanced from time  $t_0$  to time  $t^*$  for motion parallel to the geomagnetic field. These values are then used to advance the solution to the next time step  $(t_1)$  for motion perpendicular to the geomagnetic field. The splitting method is used because the equations are solved differently for the parallel and perpendicular motions. Each of these methods is described in sections 4.1 and 4.2, respectively.

# 4.1. Parallel Transport

The motion of the plasma along the geomagnetic field is described by a set of advection/diffusion equations. The advection terms are solved using an implicit donor cell method, while the diffusion terms are backward biased for stability. We illustrate this scheme in detail for the ion continuity equation:

$$\frac{\partial n}{\partial t} + b_s^2 \frac{\partial}{\partial s} \frac{nV_s}{b_s} = \mathcal{P} - n\mathcal{L} , \qquad (39)$$

where  $\mathcal{P}$  are the ion production terms and  $\mathcal{L}$  are the ion loss terms. The subscript i has been neglected for simplicity. The loss term is written as  $n\mathcal{L}$  (i.e., directly proportional to the density n) for numerical stability.

The grid is in the s direction. The physical parameters (i.e., density, velocity, and temperature) are defined at  $s_{j-1}$ ,  $s_j$ ,  $s_{j+1}$ , etc. Using this grid, the difference equation for continuity is written as

$$\frac{n_j^{t+\Delta t} - n_j^t}{\Delta t} + w \frac{(n^{t+\Delta t}V)_{j+1/2} - (n^{t+\Delta t}V)_{j-1/2}}{\Delta s_j}$$
$$= \mathcal{P} - n_j^{t+\Delta t} \mathcal{L} , \qquad (40)$$

where  $\Delta s_j = (s_{j+1} - s_{j-1})/2$  and V,  $\mathcal{P}$ , and  $\mathcal{L}$  are evaluated at time t; the subscript s on the velocity has been neglected for simplicity. Equation (40) is said to be in "flux" or "conservation" form. The density is evaluated at the upper time level  $(t + \Delta t)$  so that the difference scheme is implicit (i.e., backward biased). This method allows the Courant condition  $(\Delta t < \Delta l/V)$  to be based on the advection velocity  $(V = V_s)$  and not the sum of the advection velocity and the sound speed  $(V = V_s + C_s)$ . This technique is useful in modeling subsonic plasmas such as the ionosphere.

The donor cell method is used to compute the fluxes nV at the half steps j - 1/2 and j + 1/2. The velocity at the half step is calculated as a simple average:

$$V_l = \frac{1}{2}(V_{j-1} + V_j); \quad V_r = \frac{1}{2}(V_j + V_{j+1}) .$$
 (41)

The difference equation is written as

$$\frac{n_{j}^{t+\Delta t} - n_{j}^{t}}{\Delta t} + \frac{a_{0}n_{j-1}^{t+\Delta t} + b_{0}n_{j}^{t+\Delta t} + c_{0}n_{j+1}^{t+\Delta t}}{\Delta s_{j}}$$

$$= \mathcal{P} - n_j^{t+\Delta t} \mathcal{L} , \qquad (42)$$

where

$$a_0 = \begin{cases} -V_l & V_r > 0 & V_l > 0\\ 0 & V_r < 0 & V_l < 0\\ 0 & V_r > 0 & V_l < 0\\ -V_l & V_r < 0 & V_l > 0 \end{cases}$$
(43)

$$b_{0} = \begin{cases} V_{r} & V_{r} > 0 \ V_{l} > 0 \\ -V_{l} & V_{r} < 0 \ V_{l} < 0 \\ V_{r} - V_{l} & V_{r} > 0 \ V_{l} < 0 \\ 0 & V_{r} < 0 \ V_{l} > 0 \end{cases}$$
(44)

$$c_0 = \begin{cases} 0 & V_r > 0 & V_l > 0 \\ V_r & V_r < 0 & V_l < 0 \\ 0 & V_r > 0 & V_l < 0 \\ V_r & V_r < 0 & V_l > 0 \end{cases}$$
(45)

The main idea of the donor cell method is that the flux into or out of each cell is determined by the direction of the velocity at the cell boundary. If the velocity at the half step j-1/2 is positive (i.e.,  $V_l>0$ ), then the flux at the boundary is  $n_{j-1}V_l$  because plasma from cell j-1 is being carried into cell j. On the other hand, if the velocity at the half step j-1/2 is negative (i.e.,  $V_l<0$ ), then the flux at the boundary is  $n_jV_l$  because plasma is carried from cell j into cell j-1.

Finally, (42) is written in the form

$$A n_{j-1}^{t+\Delta t} + B n_{j}^{t+\Delta t} + C n_{j+1}^{t+\Delta t} = D , \qquad (46)$$

where

$$A = \frac{a_0}{\Delta s_j} ,$$
 
$$B = \frac{1}{\Delta t} + \frac{b_0}{\Delta s_j} + \mathcal{L} ,$$
 
$$C = \frac{c_0}{\Delta s_j} ,$$
 
$$D = n_j^t + \mathcal{P} ,$$

and A, B, C, and D are evaluated at time t. The density at the updated time step  $t + \Delta t$  is obtained by solving the tridiagonal matrix (46). The boundary condition at each end of the flux tube is  $n_i = \mathcal{P}_i/\mathcal{L}_i$ .

The momentum equation is solved in a similar manner. The advective term in the momentum equation is approximated as

$$(\mathbf{V} \cdot \nabla) V_s \simeq V_s \frac{\partial V_s}{\partial s} = \frac{1}{2} \frac{\partial V_s V_s}{\partial s} .$$

For the results presented in this paper this advective term is unimportant because the flows are subsonic, but it is included for completeness. The boundary condition at each end of the flux tube is  $V_{is} = 0$ .

The ion temperature equation is somewhat more complicated because the advection terms cannot be written entirely in conservation form and because there is a diffusion term. The ion temperature equation along the geomagnetic field is

$$\frac{\partial T_i}{\partial t} + b_s V_{is} \frac{\partial T_i}{\partial s} + \frac{2}{3} T_i b_s^2 \frac{\partial}{\partial s} \frac{V_{is}}{b_s} + \frac{2}{3} \frac{1}{n_i k} b_s^2 \frac{\partial}{\partial s} \kappa_i \frac{\partial T_i}{\partial s} 
= Q_{in} + Q_{ij} + Q_{ie} .$$
(47)

The advection terms are rewritten to obtain at least one term in conservation form

$$b_{s}V_{is}\frac{\partial T_{i}}{\partial s} + \frac{2}{3}T_{i}b_{s}^{2}\frac{\partial}{\partial s}\frac{V_{is}}{b_{s}}$$

$$= \frac{\partial b_{s}V_{is}T_{i}}{\partial s} - \frac{5}{3}V_{is}T_{i}\frac{\partial b_{s}}{\partial s} - \frac{1}{3}b_{s}T_{i}\frac{\partial V_{is}}{\partial s}. \tag{48}$$

The ion temperature equation then becomes

$$\frac{\partial T_i}{\partial t} + \frac{\partial b_s V_{is} T_i}{\partial s} - \frac{5}{3} V_{is} T_i \frac{\partial b_s}{\partial s} - \frac{1}{3} b_s T_i \frac{\partial V_{is}}{\partial s} + \frac{2}{3} \frac{1}{n_i k} b_s^2 \frac{\partial}{\partial s} \kappa_i \frac{\partial T_i}{\partial s} = Q_{in} + Q_{ij} + Q_{ie} .$$
(49)

Following similar techniques as described above, the ion temperature equation is differenced as

$$\frac{T_{j}^{t+\Delta t} - T_{j}^{t}}{\Delta t} + \frac{(b_{s}V_{s}T^{t+\Delta t})_{j-1/2} - (b_{s}V_{s}T^{t+\Delta t})_{j-1/2}}{\Delta s_{j}}$$

$$-\frac{5}{3}V_{s}T_{j}^{t+\Delta t}\frac{b_{s,j+1} - b_{s,j-1}}{\Delta s'_{j}} - \frac{1}{3}b_{s}T_{j}^{t+\Delta t}\frac{V_{s,j+1} - V_{s,j-1}}{\Delta s'_{j}}$$

$$-\frac{2}{3}\frac{b_{s}^{2}}{nk} \left[\frac{\kappa_{j+1/2}}{\Delta s'_{j}}\frac{T_{j+1}^{t+\Delta t} - T_{j}^{t+\Delta t}}{ds_{j}}\right]$$

$$-\frac{\kappa_{j-1/2}}{\Delta s'_{j}}\frac{T_{j}^{t+\Delta t} - T_{j-1}^{t+\Delta t}}{ds_{j-1}}$$

$$= [Q_{in} + Q_{ij} + Q_{ie}] - [Q'_{in} + Q'_{ij} + Q'_{ie}]T_{j}^{t+\Delta t}, (50)$$

where

$$\begin{split} \mathcal{Q}_{in} &= \sum_{q} \frac{2m_{i}m_{q}}{\left(m_{i} + m_{q}\right)^{2}} \nu_{iq} \left[ kT_{q} + \frac{1}{3}m_{q} \left| V_{qs} - V_{is} \right|^{2} \right] \;, \\ \mathcal{Q}_{ij} &= \sum_{j} \frac{2.2 \times 10^{-4} n_{j}}{A_{i}A_{j} \left( T_{i}/A_{i} + T_{j}/A_{j} \right)^{3/2}} T_{j} \;, \\ \mathcal{Q}_{ie} &= \frac{7.7 \times 10^{-6} n_{e}}{A_{i}T_{e}^{3/2}} T_{e} \;, \\ \mathcal{Q}'_{in} &= \sum_{q} \frac{2m_{i}m_{q}}{\left(m_{i} + m_{q}\right)^{2}} \nu_{iq} k \;, \\ \mathcal{Q}'_{ij} &= \sum_{j} \frac{2.2 \times 10^{-4} n_{j}}{A_{i}A_{j} \left( T_{i}/A_{i} + T_{j}/A_{j} \right)^{3/2}} \;, \\ \mathcal{Q}'_{ie} &= \frac{7.7 \times 10^{-6} n_{e}}{A \cdot T^{3/2}} \;. \end{split}$$

The velocity at the half-step is calculated as a simple average but now includes the magnetic flux factor  $b_s$ 

$$V_l = \frac{1}{2}(b_{s,j-1}V_{j-1} + b_{s,j}V_j)$$
 (51)

$$V_r = \frac{1}{2} (b_{s,j} V_j + b_{s,j+1} V_{j+1}) . {(52)}$$

(53)

Finally, (50) is written in the form

$$AT_{i-1}^{t+\Delta t} + BT_{i}^{t+\Delta t} + CT_{i+1}^{t+\Delta t} = D$$
, (54)

where

$$A = \frac{a_0}{\Delta s_j} - \frac{2}{3} \frac{b_s^2}{nk} \frac{\kappa_{j-1/2}}{\Delta s_j'} \frac{1}{ds_{j-1}} ,$$

$$B = \frac{1}{\Delta t} + \frac{b_0}{\Delta s_j} - \frac{5}{3} V_s \frac{b_{s,j+1} - b_{s,j-1}}{\Delta s_j'}$$

$$- \frac{1}{3} b_s \frac{V_{s,j+1} - V_{s,j-1}}{\Delta s_j'}$$

$$+ \frac{2}{3} \frac{b_s^2}{nk} \left[ \frac{\kappa_{j+1/2}}{\Delta s_j'} \frac{1}{ds_j} + \frac{\kappa_{j-1/2}}{\Delta s_j'} \frac{1}{ds_{j-1}} \right]$$

$$+ [\mathcal{Q}'_{in} + \mathcal{Q}'_{ij} + \mathcal{Q}'_{ie}] ,$$

$$C = \frac{c_0}{\Delta s_j} - \frac{2}{3} \frac{b_s^2}{nk} \frac{\kappa_{j+1/2}}{\Delta s_j'} \frac{1}{ds_j} ,$$

$$D = T_{i,j}^t + [\mathcal{Q}_{in} + \mathcal{Q}_{ij} + \mathcal{Q}_{ie}] .$$

In the above, A, B, C, and D are evaluated at time t, and  $a_0$ ,  $b_0$ , and  $c_0$  are defined as in (43) – (45) but with the appropriate velocities (i.e., equation (53)). The boundary condition at each end of the flux tube is  $T_i = T_{-}$ 

The electron temperature equation (31) is solved the same way as the ion temperature equation but with no advection terms (i.e.,  $V_{es} = 0$ ) and the appropriate heating terms (equations (32) – (34)).

#### 4.2. Perpendicular Transport

The perpendicular transport scheme involves two stages. First, we calculate the new value of p that the flux tube drifts to because of the  $\mathbf{E} \times \mathbf{B}$  drift. Second, we update the density and temperature of each species at the new value of p based on conservation of particle number and magnetic flux. We assume that there is no flux transport along  $\mathbf{B}$  during this stage.

Consider two flux tubes at times t and  $t + \Delta t$ , schematically shown in Figure 2. The flux tube at time t has  $p = p_1$ ; it rises to  $p = p_2$  at time  $t + \Delta t$  because of an  $\mathbf{E} \times \mathbf{B}$  drift  $V_E$ . At time t,  $p_1$  and the s coordinates

are known, as well as the physical variables along  $p_1$ . The value of  $p_2$  is determined as follows. From (2) we know that

$$\Delta p = p_2 - p_1 = (r_{e2} - r_{e1})/R_E = \Delta r_e/R_E$$
, (55)

along s=0 (i.e.,  $\theta_e=\pi/2$ ). Since  $\Delta r_e=V_E\Delta t$ , we find that

$$p_2 = p_1 + V_E \Delta t / R_E. \tag{56}$$

Finally, we note that the magnitude of the geomagnetic field is also known along  $p = p_2$ .

We now calculate the updated density and temperature along  $p=p_2$ . The volume of flux tube cell at time t is  $\mathcal{V}_1=\pi\rho_1^2dl_1$ , and at time  $t+\Delta t$  it is  $\mathcal{V}_2=\pi\rho_2^2dl_2$ , where r is the radius of the flux tube and dl is the length of the flux tube along the geomagnetic field. This is illustrated in Figure 3. The total number of particles in each flux tube cell is conserved so that

$$n_1 \rho_1^2 dl_1 = n_2 \rho_2^2 dl_2. (57)$$

The magnetic flux through each flux tube is also conserved so that

$$\rho_1^2 B_1 = \rho_2^2 B_2. \tag{58}$$

Combining (57) and (58), we find that the updated density  $n_2$  at time  $t + \Delta t$  is given by

$$n_2 = \frac{dl_1}{dl_2} \frac{B_2}{B_1} n_1 \ . \tag{59}$$

The updated temperature is given by

$$T_2 = \frac{dl_1}{dl_2} \frac{B_2}{B_1} T_1 \ , \tag{60}$$

where it is assumed that  $\gamma = 1$  in the equation of state.

A key feature of this scheme is that the s coordinate is the same for all values of p that the flux tube may have. We do this in SAMI2 as follows. For a given initial flux tube we calculate the maximum p value that the flux tube can have during a 24-hour period using a prescribed electric field model. We then define the s coordinates for this flux tube using the method described in section 2. We typically use 201 grid points with the first and last grid points at an altitude of 90 km. When the flux tube is not at its maximum altitude, a number of the grid points at each end of the flux tube will fall below 90 km as shown in Figure 4; the solid lines in Figure 4 show the mapping of s values from the flux tube at a low altitude to a higher altitude. Only grid points above 90 km are used in the numerical algorithms to update the variables. The boundary conditions at each end of the flux tube (i.e., altitude of  $\sim 90 \text{ km}$ ) are  $n_i = \mathcal{P}_i/\mathcal{L}_i$ ,  $T_i = T_n$ , and  $V_{is} = 0$ .

Finally, one can incorporate a drift in the  $\phi$  direction in a straightforward manner. One calculates the new longitudinal coordinate that the flux tube drifts to using

$$\phi_2 = \phi_1 + \sin^{-1}(\Delta x/r_{e2}) , \qquad (61)$$

where  $\Delta x = V_{\phi} \Delta t$ .

# 5. Photodeposition Model

# 5.1. Daytime Model

The daytime photoionization rate of neutrals is given by

$$q_l(z) = n_l(z) \sum_{\lambda} \sigma_l^{(i)}(\lambda) \phi_{\infty}(\lambda) \times$$

$$\exp\left[-\sum_{m} \sigma_{m}^{(a)}(\lambda) \int_{z}^{\infty} n_{m}(s) ds\right] , \qquad (62)$$

where  $q_l(z)$  is the photoionization rate of the lth neutral species at altitude z,  $n_{l,m}(z)$  is the neutral density of l, m at  $z, \phi_{\infty}(\lambda)$  is the incident solar flux at wavelength  $\lambda$ ,  $\sigma_l^{(i)}(\lambda)$  is the photoionization cross-section for l at  $\lambda$ , and  $\sigma_m^{(a)}(\lambda)$  is the photoabsorption cross-section for m at  $\lambda$ 

The integral  $\int_{z}^{\infty} n_{m}(s) ds$  is evaluated along a ray from the Sun to z. The ray is described by the zenith angle  $\chi$ , which is defined by the angle that the ray makes with the vertical through z. In the original SAMI this integral is actually calculated numerically. Here we use a simplified form based on an exponential neutral atmosphere [Smith and Smith, 1972]. This method is not highly accurate for grazing solar zenith angles.

Chapman [1931] writes the columnar content as

$$I_p = n_p H \operatorname{ch}(X_p, \chi_p) \tag{63}$$

where

$$I_p = \int_{s_n}^{\infty} n_n(s) \, ds \;, \tag{64}$$

and it is assumed that

$$n_n = n_{n0} \exp(-h/H) , \qquad (65)$$

where h is the altitude  $(h = h_p \text{ at point } P)$ , H is the scale height  $(H = kT_n/m_ng)$ ,  $n_n$  is the density of neutral species  $(n_n = n_{n0} \text{ at } h = 0)$ ,  $n_p$  is the density  $n_n$  at point P,  $r_p$  is the radial distance to point P, s is the distance coordinate measured positively in the solar direction  $(s = s_p \text{ at point } P)$ , and  $X_p = r_p/H$ . The geometry for the photodeposition model is shown in Figure 5. The form of  $\operatorname{ch}(X_p, \chi_p)$  depends on whether or not the solar zenith angle  $\chi_p$  at point P is less than

or greater than  $90^{\circ}$ . The current version of SAMI2 only considers  $\chi_p < 90^{\circ}$ . For this case we use

$$\operatorname{ch}(X_p, \chi_p) \simeq \left[ (\pi/2) X_p \right]^{1/2} \exp(y^2) (1 - \operatorname{erfy}) ,$$
 (66)

where

$$y = \left(X_p/2\right)^{1/2} \left|\cos \chi_p\right| \, .$$

and erf is the error function. Thus we find that

$$I_p \simeq n_p H_p \left[ (\pi/2) X_p \right]^{1/2} \exp(y^2) \cdot (1 - \text{erfy})$$
 (67)

The solar EUV flux  $\phi_{\infty}(\lambda)$  is determined from the EUVAC model developed by *Richards et al.* [1994]. The flux is determined in 37 wavelength bins (denoted by the subscript i) using

$$\phi_i = F74113_i \left[ 1 + A_i (P - 80) \right] , \qquad (68)$$

where  $F74113_i$  is the modified reference flux,  $A_i$  is the scaling factor, and P = (F10.7A + F10.7)/2 (F10.7A is an 81-day average of the daily F10.7 index). The values used in SAMI2 are listed in Table 1. We consider the photoionization of four species: He, O, N<sub>2</sub>, and O<sub>2</sub>. The photoionization cross sections used are given in Table 2.

#### 5.2. Nighttime Model

Nighttime photoionization of neutrals is given by

$$q_l(z) = n_l(z) \sum_{\lambda} \sigma_l^{(i)}(\lambda) \phi(z, \lambda) , \qquad (69)$$

where  $q_l(z)$  is the photoionization rate of the lth neutral species at altitude z,  $n_l(z)$  is the density of l at z,  $\phi(z, \lambda)$  is the solar flux at the altitude z and the wavelength  $\lambda$ , and  $\sigma_l^{(i)}(\lambda)$  is the photoionization cross-section for l at  $\lambda$ 

The EUV lines and photoionization cross sections considered are given in Table 3, where the photoionization cross sections are in units  $10^{-18}$  cm<sup>2</sup>; they are obtained from the original SAMI code data files [Oran et al., 1974]. The EUV flux is obtained from Strobel et al. [1974]. The EUV flux is given for each line as a function of altitude and solar zenith angle over limited ranges. The fluxes are approximated with analytical functions, which we use to interpolate fluxes for all altitudes and zenith angles in the range  $90^{\circ} < \theta < 180^{\circ}$ .

# 6. Chemistry Model

The chemical reactions and reaction rates used in SAMI2 are given in Table 4; the recombination loss

rates are given in Table 5. These provide source  $(\mathcal{P}_i)$  and loss  $(\mathcal{L}_i)$  terms in the continuity equation.

The chemical reactions are implemented in the continuity equation as follows. Consider the general chemical reaction

$$X^+ + Y \to X + Y^+ \tag{70}$$

that occurs at a rate  $\kappa_{XY}$ . This reaction represents a loss of X ions and a source of Y ions. The loss and production rates of X and Y are the same:

$$\mathcal{L}_X = \kappa_{XY} n_i(X) n_n(Y) = \mathcal{P}_Y , \qquad (71)$$

where  $n_i(X)$  is the density of X ions and  $n_n(Y)$  is the density of Y neutrals. However, we have factored out the ion density in defining the loss term  $L_i$  in the continuity equation. Thus the chemical loss term used in SAMI2 is simply

$$\mathcal{L}_X = \kappa_{XY} n_n(Y) \ . \tag{72}$$

The radiative recombination loss of ions is represented by the reaction

$$X^+ + e \to X \tag{73}$$

at a rate  $\kappa_{Xe}$ . The loss of  $X^+$  is given by

$$\mathcal{L}_X = \kappa_{Xe} n_i(X) n_e , \qquad (74)$$

where  $n_e$  is the electron density. The recombination loss rate used in SAMI2 is then

$$\mathcal{L}_X = \kappa_{Xe} n_e \ . \tag{75}$$

#### 7. Results

The simulation results presented use a mesh of 201 grid points for each flux tube, and 120 flux tubes are modeled. The time step used is  $\Delta t = 1-12$  s, depending upon the altitude of the flux tube. The lower boundary of each flux tube is 90 km. The simulations are run for 24 hours starting shortly after sunrise (t=0700 LT). A model  $\mathbf{E} \times \mathbf{B}$  drift is used that is proportional to  $\sin[(t-7)/24]$ , where t is the local time in hours. In general, the ionospheric flux tubes rise during the day from 0700 to 1900 LT, and they fall during the night from 1900 to 0700 LT. In addition to the temporal dependence, we also include an altitude dependence in the  $\mathbf{E} \times \mathbf{B}$  drift velocity. Specifically, we take

$$V_E(z) = \begin{cases} 0 & z < 150 \text{ km} \\ V_0(z - 150)/150 & 150 \text{ km} < z < 300 \text{ km} \\ V_0 & 400 \text{ km} < z \end{cases},$$
(76)

where z is the altitude of the flux tube at s=0, i.e., at the magnetic equator. In the simulations presented here we use  $V_0=20~{\rm m\,s^{-1}}$ . This simplistic model is used so that the flux tubes do not drift out of the lower altitude range; it will be replaced in the near future with an empirical electric field model (L. Scherliess, private communication, 2000).

The first set of simulation results use the following parameters: year 1994, day 274 (October 1), longitude 283.13° E (Jicamarca), Ap 4, F10.7 75, and F10.7A 75. In Figure 6 we plot the density of atomic ions (H<sup>+</sup>, He<sup>+</sup>, and O<sup>+</sup>) and the electron density as a function of altitude over Jicamarca ( $\theta_g \simeq -12^{\circ}$ ) at times 1100, 1500, 1900, and 0200 LT. During the day a relatively broad F layer forms between 200 and 600 km with a peak density  $\sim 10^6$  cm<sup>-3</sup>. At night the ionosphere decays with a strong erosion on the bottomside, below 200 km. The peak density drops to  $\sim 10^5$  cm<sup>-3</sup>. The transition altitude for  $O^+/H^+$  ranges from  $\sim 550$  km in the early morning to  $\sim 800$  km in the late afternoon. Finally, He<sup>+</sup> is a minor ion that is at most  $\sim 10\%$  of the electron density but is typically  $\leq 4\%$  of the electron density.

In Figure 7 we plot  $N_mF_2$  and  $h_mF_2$  as a function of time over Jicamarca. The peak electron density is fairly constant during the day at  $\sim 10^6$  cm<sup>-3</sup>. At night it monotonically decreases, with the largest decrease occurring in the early morning hours. Before sunrise the peak electron density is a few  $\times 10^4$  cm<sup>-3</sup>. The  $h_mF_2$  shows more variability than the  $N_mF_2$  does. The  $h_mF_2$  rises from 240 km at sunrise to almost 400 km by late morning and remains near this altitude until 2000 LT. It then decreases to  $\sim 250$  km by early morning (0200 LT); subsequently it rises slightly to  $\sim 270$  km just before sunrise.

In Plate 1 we show a colored contour plot of the logarithm of the electron density as a function of geographic latitude versus altitude at times 1100, 1500, 1900, and 0200 LT. During the day a broad region of intense ionization develops between 200 and 500 km. By midafternoon (t=1500 LT), ionization crests are forming away from the geomagnetic equator ( $\theta_g \simeq -9^{\circ}$ ). At time 1900 LT the Appleton anomaly is evident with the maximum electron density occurring  $\pm 15^{\circ}$  from the geomagnetic equator. Also, the electron density has become very low at altitudes below 200 km. Finally, at time 0200 LT the ionosphere has decayed substantially. The electron density has fallen below  $10^3$  cm<sup>-3</sup> in the altitude range 110 – 200 km.

In Figure 8 we plot the ion O<sup>+</sup> temperature and the electron temperature as a function of altitude over Ji-

camarca at times t=0800 and 1100 LT. In the early morning, photoelectron heating rapidly heats the electrons at  $\sim 200$  km and above 300 km. At higher altitudes the electrons are hotter than the ions. However, by late morning (t=1100 LT) the electrons and ions have cooled off substantially at higher altitudes and have equilibrated. However, at lower altitudes (around 200 km) the electrons remain hotter than the ions.

The second set of simulation results use the following parameters: year 1994, day 274 (October 1), longitude 293.4° E (Arecibo), Ap 4, F10.7 75, and F10.7A 75. These are the same parameters as those used in the first set of simulation results, but they are at a different longitude. In Figure 9 we plot density as a function of altitude over Arecibo ( $\theta_g \sim 18.5^{\circ}$ ) at four times: 1100, 1500, 1900 and 0200 LT. We show the densities of the atomic ions H<sup>+</sup>, He<sup>+</sup>, and O<sup>+</sup>, as well as the electron density. The characteristics of these profiles are similar to those shown in Figure 6 over Jicamarca with some differences. During the day the topside F region does not become as broad as it does over Jicamarca. The transition altitude for  $O^+$ -H<sup>+</sup> is  $\sim 500$  km at 0200 LT and rises to  $\sim 1000$  km in the late morning. This result is consistent with the modeling results presented by MacPherson et al. [1998]. Helium ions are a minor species; the He<sup>+</sup> density peaks in the range 450 -750 km depending on the time of day. The ratio of the  $\mathrm{He^{+}}$  density to the electron density is < 8%. The altitude variation of the He<sup>+</sup> density is different than that over Jicamarca; the peak He<sup>+</sup> densities occur at lower altitudes over Arecibo than Jicamarca.

In Figure 10 we plot  $N_m F_2$  and  $h_m F_2$  as a function of time over Arecibo. The behavior of  $N_m F_2$  and  $h_m F_2$ over Arecibo is somewhat different than that over Jicamarca. The  $N_m F_2$  increases monotonically from sunrise until  $\sim 1600$  LT and then drops by an order of magnitude within  $\sim 4$  hours; it then plateaus until the early morning ( $\sim 0200$  LT) before decreasing slightly until sunrise. The  $h_m F_2$  exhibits similar behavior until  $\sim$ 1900 LT. However, during the evening from 1900 LT until midnight, the  $h_m F_2$  rises from 260 to 290 km; it then decreases somewhat to  $\sim 270$  km before sunrise. These results are consistent with those obtained from SUPIM; for example, they compare favorably to the plus curves in Figures 3 and 5 of it MacPherson et al. [1998]; the differences in the  $h_m F_2$  curves are probably due to the different electric field models used.

In Figure 11 we plot the density of the molecular ions as a function of altitude over Arecibo for the same times as those in Figure 9. The molecular ions shown are  $N_2^+$ ,  $NO^+$ , and  $O_2^+$ . The  $N_2^+$  ions only have an

appreciable density ( $\lesssim 10^3~{\rm cm}^{-3}$ ) during the daytime. In general, the densities of NO<sup>+</sup> and O<sub>2</sub><sup>+</sup> are comparable with O<sub>2</sub><sup>+</sup>, dominating in the E region in the daytime. However, observations indicate that NO<sup>+</sup> is the major ion in the E region although most models find O<sub>2</sub><sup>+</sup> to be the dominant ion [Richards and Torr, 1996]. We find that NO<sup>+</sup> is the dominant ion during the nighttime.

In Figure 12 we plot the ion  ${\rm O}^+$  temperature and the electron temperature as a function of altitude over Arecibo at times t=800 and 1100 LT. In the early morning, photoelectron heating rapidly heats the electrons at above 200 km. At higher altitudes the electrons are somewhat hotter than the ions. By late morning the ion and electron temperatures have equilibrated at altitudes above 600 km. Overall, these results are similar to those over Jicamarca.

# 8. Summary

In this paper we have described the new low-latitude ionospheric model that has been developed at the Naval Research Laboratory: Sami2 is Another Model of the Ionosphere (SAMI2). SAMI2 treats the dynamic plasma and chemical evolution of seven ion species (H<sup>+</sup>,  $He^+$ ,  $N^+$ ,  $O^+$ ,  $N_2^+$ ,  $NO^+$ , and  $O_2^+$ ) in the altitude range  $\sim 100~\mathrm{km}$  to several thousand kilometers. The ion continuity and momentum equations are solved for all seven species; the temperature equation is solved for H<sup>+</sup>, He<sup>+</sup>, O<sup>+</sup>, and the electrons. SAMI2 models the plasma along the Earth's dipole field from hemisphere to hemisphere, includes the  $\mathbf{E} \times \mathbf{B}$  drift of a flux tube (both in altitude and in longitude), and includes ion inertia in the ion momentum equation for motion along the dipole field line. The final point is relevant for plasma dynamics at very high altitudes, where ion inertia becomes important. The neutral species are modeled using the empirical models MSIS86 and HWM93.

We have also implemented a parallelized version SAMI2 using the Message Passing Interface (MPI) methods. This method leads to a greatly reduced runtime based on the number of processors available. For example, simulations of a 24 hour period using 120 flux tubes with 201 grid points can be completed within 1 hour using 10 Pentium 550-MHz processors.

Recently, several new results have been reported based on simulation studies using SAMI2. First, we have discovered the possibility of the formation of an electron hole (depletion) in the topside equatorial ionosphere [ $Huba\ et\ al.$ , 2000a]. The reduction in the electron density occurs in the altitude range 1500 - 2500 km at geomagnetic equatorial latitudes. The hole is produced by transhemispheric  $O^+$  flows that collision-

ally couple to H<sup>+</sup> and transport it to lower altitudes, thereby reducing the electron density at high altitudes. The transhemispheric O<sup>+</sup> flows are caused by an interhemispheric pressure anisotropy that can be generated by the neutral wind, primarily during solstice conditions. The formation of the electron hole has a seasonal and longitudinal dependence. Second, we have also found that ion sound waves may be generated in the topside low-latitude ionosphere [Huba et al., 2000b]. The waves are excited at sunrise and sunset; they persist  $\sim 1-3$  hours at altitudes above  $\sim 1000$  km with periods on the order of tens of minutes. The waves are the result of the rapid heating and cooling of the lower ionosphere that occurs at sunrise and sunset. At sunrise, photoelectron heating produces strong upward plasma flows along the geomagnetic field. These flows lead to a local compression and heating of the plasma at the apex of the field line which, in turn, generates ion sound waves. At sunset the waves are produced by a rapid cooling of the plasma. This result is obtained only because ion inertia is included in SAMI2.

A number of improvements are planned for SAMI2. First, a simple analytical model is used to model the zonal electric field drift; this is done for simplicity and to assure that we are able to recover the Appleton anomaly. However, most low-latitude models use an empirical electric field model[Fejer, 1991,1997] that is more realistic than the simple model used in SAMI2. Recently, Scherliess and Fejer [1999] developed a global equatorial F region vertical drift model based on radar and satellite observations. In the near future this model will be available as a FORTRAN code, similar to MSIS and HWM, and will be incorporated into SAMI2 (L. Scherliess, private communication, 2000). Second, we are also developing a self-consistent formulation of the zonal electric field using MSIS and HWM. This involves solving a potential equation for the electrostatic potential in two dimensions. The thermospheric winds provide the source of energy to maintain the zonal electric field; an excellent overview of the dynamo process is presented by Kelley [1989]. Third, the photoelectron heating model will be improved. A more physically based photoelectron heating model has been developed by Richards and Torr [1996] based on the two-stream model of Nagy and Banks [1970]. We intend to incorporate this model into SAMI2 in the near future. Finally, there are a number of additional heating terms that will be included in the energy equations such as elastic collisions between O<sup>+</sup>, H<sup>+</sup>, and He<sup>+</sup> and neutrals, and electron heating associated with fine structure [Bailey and Balan, 1996].

# Appendix A: Coordinate Transformations

The notation for the subscripts is as follows: g, geographic spherical, t, tilted spherical, e, eccentric spherical, and d, dipole. This analysis is similar to that presented by Bailey and Balan [1996] and Millward et al. [1996].

#### A1. Spherical Geographic to Spherical Tilted

The geometry associated with th transformation from spherical geographic to spherical tilted is shown in Figure 13. The angles that define the tilted dipole system are  $\theta_n=11^0$  and  $\phi_n=289^0$  [Millward et al., 1996]. A general point P is defined by the coordinates  $(\theta_g,\phi_g,r_g)$  in the spherical geographic system where  $N_g$  is the geographic North Pole. We need to determine the coordinates of P in the spherical tilted system  $(\theta_t,\phi_t,r_t)$  where  $N_t$  is the magnetic North Pole in the tilted dipole system. It is clear that  $r_g=r_t$  since the two systems have the same origin. The angles  $\theta_t$  and  $\phi_t$  can be determined using the formulas for spherical triangles. Referring to Figure 14a we note that

$$\cos a = \cos b \cos c + \sin b \sin c \cos A . \tag{77}$$

By making the corresponding identifications with the angles defined in Figure 13 (shown in Figure 14b) and using (77), we find the following relationships:

$$\sin \theta_t = \sin \theta_q \sin \theta_n + \cos \theta_q \cos \theta_n \cos (\phi_q - \phi_n) , (78)$$

$$\sin \theta_a = \sin \theta_t \sin \theta_n - \cos \theta_t \cos \theta_n \cos \phi_t , \qquad (79)$$

$$\cos a = \cos b \cos c + \sin b \sin c \cos A . \tag{80}$$

In deriving (78) – (80) we used the relations  $\sin(\pi/2 - \theta) = \sin \theta$ ,  $\cos(\pi/2 - \theta) = \cos \theta$ , and  $\cos(\pi - \theta) = -\cos \theta$ . The transformation equations from  $(\theta_g, \phi_g)$  to  $(\theta_t, \phi_t)$  are

$$\theta_t = \sin^{-1} \left[ \sin \theta_g \sin \theta_n + \cos \theta_g \cos \theta_n \cos (\phi_g - \phi_n) \right] ,$$
(81)

$$\phi_t = \begin{cases} \phi_{t0} & \phi_n - \pi > \phi_g > \phi_n \\ 2\pi - \phi_{t0} & \phi_n - \pi < \phi_g < \phi_n \end{cases}, \tag{82}$$

$$\phi_{t0} = \cos^{-1} \left[ \frac{\sin \theta_t \sin \theta_n - \sin \theta_g}{\cos \theta_t \cos \theta_n} \right] . \tag{83}$$

The reverse transformation equations from  $(\theta_t, \phi_t)$  to  $(\theta_q, \phi_q)$  are

$$\theta_g = \sin^{-1} \left[ \sin \theta_t \sin \theta_n - \cos \theta_t \cos \theta_n \cos \phi_t \right] , \quad (84)$$

$$\phi_g = \begin{cases} \phi_n + \phi_{g0} - 2\pi & \phi_t < \pi \\ \phi_n - \phi_{g0} & \phi_t > \pi \end{cases}$$
 (85)

$$\phi_{g0} = \cos^{-1} \left[ \frac{\sin \theta_t - \sin \theta_g \sin \theta_n}{\cos \theta_g \cos \theta_n} \right] , \qquad (86)$$

and  $\phi_g = \phi_g + 2\pi$  if  $\phi_g < 0$ .

### A2. Spherical Tilted to Spherical Eccentric

The transformation from the spherical tilted system to the spherical eccentric system only involves a translation. This is shown in Figure 15 where the displacement of the two systems is defined by  $(x_{0t}, y_{0t}, z_{0t})$ . The Cartesian coordinates of point P in the tilted system are

$$x_t = r_t \cos \theta_t \cos \phi_t , \qquad (87)$$

$$y_t = r_t \cos \theta_t \sin \phi_t \,\,, \tag{88}$$

$$z_t = r_t \sin \theta_t , \qquad (89)$$

and in the eccentric system they are

$$x_e = x_t - x_{0t} = r_t \cos \theta_t \cos \phi_t - x_{0t}$$
, (90)

$$y_e = y_t - y_{0t} = r_t \cos \theta_t \sin \phi_t - y_{0t} , \qquad (91)$$

$$z_e = z_t - z_{0t} = r_t \sin \theta_t - z_{0t} . {92}$$

The transformation from spherical tilted to spherical eccentric coordinates is given by

$$r_e = (x_e^2 + y_e^2 + z_e^2)^{1/2} ,$$
 (93)

$$\theta_e = \tan^{-1}(z_e/\rho_e) , \qquad (94)$$

$$\phi_e = \begin{cases} \phi_{e0} & \phi_{e0} > 0\\ \phi_{e0} + 2\pi & \phi_{e0} < 0 \end{cases}, \tag{95}$$

where  $\rho_e = \left(x_e^2 + y_e^2\right)^{1/2}$ ,  $\phi_{e0} = \tan^{-1}\left(y_e/x_e\right)$ , and one substitutes (90) - (92) into (93) - (95).

The reverse transformation from spherical eccentric to spherical tilted follows in a similar manner. The displacement of the two systems is defined by  $(x_{0g}, y_{0g}, z_{0g})$ , which corresponds to the origin of the spherical system, i.e.,  $(-x_{0t}, -y_{0t}, -z_{0t})$ . The Cartesian coordinates of point P in the eccentric system are

$$x_e = r_e \cos \theta_e \cos \phi_e , \qquad (96)$$

$$y_e = r_e \cos \theta_e \sin \phi_e \,\,, \tag{97}$$

$$z_e = r_e \sin \theta_e , \qquad (98)$$

and in the tilted system they are

$$x_t = x_e + x_{0t} = r_e \cos \theta_e \cos \phi_e + x_{0t}$$
, (99)

$$y_t = y_e + y_{0t} = r_e \cos \theta_e \sin \phi_e + y_{0t}$$
, (100)

$$z_t = z_e + z_{0t} = r_e \sin \theta_e + z_{0t} . {101}$$

The transformation from spherical eccentric to spherical tilted is given by

$$r_t = (x_t^2 + y_t^2 + z_t^2)^{1/2} ,$$
 (102)

$$\theta_t = \tan^{-1}(z_t/\rho_t) , \qquad (103)$$

$$\phi_t = \begin{cases} \phi_{t0} & \phi_{t0} > 0\\ \phi_{t0} + 2\pi & \phi_{t0} < 0 \end{cases}$$
 (104)

where  $\rho_t = (x_t^2 + y_t^2)^{1/2}$ ,  $\phi_{t0} = \tan^{-1}(y_t/x_t)$ , and one substitutes (99) - (101) into (102) - (104).

# Appendix B: Dipole Curvilinear Coefficients

The curvilinear factors h are given by [Orens et al., 1979]

$$h_q = \frac{r^3}{R_E^2} \frac{1}{\left[1 + 3\cos^2\theta\right]^{1/2}} , \qquad (105)$$

$$h_p = \frac{R_E \sin^3 \theta}{\left[1 + 3\cos^2 \theta\right]^{1/2}} , \qquad (106)$$

$$h_{\phi} = r \sin \theta \ . \tag{107}$$

The divergence is defined as

$$\nabla \cdot \mathbf{A} = \frac{1}{h_p h_q h_\phi} \left[ \frac{\partial}{\partial p} A_p h_q h_\phi + \frac{\partial}{\partial q} A_q h_p h_\phi + \frac{\partial}{\partial \phi} A_\phi h_p h_q \right] , \qquad (108)$$

and the gradient is defined as

$$\nabla f = \frac{1}{h_p} \frac{\partial f}{\partial p} \mathbf{e}_p + \frac{1}{h_q} \frac{\partial f}{\partial q} \mathbf{e}_q + \frac{1}{h_\phi} \frac{\partial f}{\partial \phi} \mathbf{e}_\phi . \tag{109}$$

Substituting in we get

$$\nabla f = \frac{\left[1 + 3\cos^2\theta\right]^{1/2}}{R_E \sin^3\theta} \frac{\partial f}{\partial p} \mathbf{e}_p +$$

$$\frac{R_E^2 \left[1 + 3\cos^2\theta\right]^{1/2}}{r^3} \frac{\partial f}{\partial q} \mathbf{e}_q + \frac{1}{r\sin\theta} \frac{\partial f}{\partial \phi} \mathbf{e}_\phi , \qquad (110)$$

$$h_p h_q h_\phi = \frac{r^4 \sin^4 \theta}{R_E \left[ 1 + 3 \cos^2 \theta \right]} , \qquad (111)$$

$$\frac{1}{h_p h_q h_\phi} = \frac{R_E \left[ 1 + 3\cos^2\theta \right]}{r^4 \sin^4\theta} \tag{112}$$

$$\nabla \cdot \mathbf{A} = \frac{\left[1 + 3\cos^2\theta\right]}{R_E r^4 \sin^4\theta} \frac{\partial}{\partial p} \frac{r^4 \sin\theta A_p}{\left[1 + 3\cos^2\theta\right]^{1/2}}$$

$$+\frac{R_E^2\left[1+3\cos^2\theta\right]}{r^4\sin^4\theta}\frac{\partial}{\partial q}\frac{r\sin^4\theta A_q}{\left[1+3\cos^2\theta\right]^{1/2}}+\frac{1}{r\sin\theta}\frac{\partial A_\phi}{\partial \phi}$$
(113)

# Appendix C: Collisional Coefficients

The ion thermal conductivity  $K_i$  is given by [Banks and Kockarts, 1973]

$$K_i = 1.24 \times 10^4 \frac{NT_i^{5/2}}{n_e} ,$$
 (114)

$$N = n(O^{+}) + 2n(He^{+}) + 4n(H^{+}), \qquad (115)$$

where  $K_i$  is in  $\operatorname{ev} \operatorname{cm}^{-1} \operatorname{s}^{-1} \operatorname{K}^{-1}$ .

The electron thermal conductivity  $K_e$  is given by [Banks and Kockarts, 1973; Bailey and Balan, 1996]

$$K_e = \frac{7.7 \times 10^5 T_e^{5/2}}{1 + 3.22 \times 10^4 T_e^2 N_q / n_e} , \qquad (116)$$

where

$$N_q = n(O)q(O) + n(N_2)q(N_2) + n(O_2)q(O_2) \; , \ \, (117)$$

$$q(0) = 1.1 \times 10^{-16} \left( 1 + 5.7 \times 10^{-4} T_e \right) ,$$
 (118)

$$q(N_2) = 2.82 \times 10^{-17} T_e^{1/2} \left( 1 - 1.21 \times 10^{-4} T_e \right) , (119)$$

$$q(O_2) = 2.2 \times 10^{-16} \left( 1 + 3.6 \times 10^{-2} T_e^{1/2} \right) ,$$
 (120)

and  $K_e$  is in ev cm<sup>-1</sup>s<sup>-1</sup>K<sup>-1</sup>.

The ion-ion collision frequency is given by [Huba, 1998]

$$\nu_{ij} = 6.8 \times 10^{-8} n_j Z_i^2 Z_j^2 \lambda_{ij} \frac{A_j^{1/2}}{A_i} \left( 1 + \frac{A_j}{A_i} \right)^{-1/2} \frac{1}{T^{3/2}} ,$$
(121)

where  $A_i$  and  $A_j$  are the atomic ion numbers, T is in eV, Z is the charge, and

$$\lambda_{ij} = 23 - \ln \left[ \frac{Z_i Z_j (A_i + A_j)}{A_i T_j + A_j T_i} \left( \frac{n_i Z_i^2}{T_i} + \frac{n_j Z_j^2}{T_j} \right)^{1/2} \right].$$
(122)

We convert from ev to kelvins using

$$T_{eV} = 8.6174 \times 10^{-5} T_K$$
, (123)

and we assume singly ionized ions, i.e.,  $Z_i = Z_j = 1$ . We find that

$$\nu_{ij} = 9.2 \times 10^{-2} n_j \lambda_{ij} \frac{A_j^{1/2}}{A_i} \left( 1 + \frac{A_j}{A_i} \right)^{-1/2} \frac{1}{T^{3/2}} ,$$

$$\lambda_{ij} = 23 - 1.25 \times 10^6 \ln \left[ \frac{A_i + A_j}{A_i T_j + A_j T_i} \left( \frac{n_i}{T_i} + \frac{n_j}{T_j} \right)^{1/2} \right] .$$
(125)

The ion-neutral collision frequency is given by [Banks and Kockarts, 1973]

$$\nu_{in} = \frac{m_n}{m_i + m_n} \bar{\nu}_{in} , \qquad (126)$$

$$\bar{\nu}_{in} = 2.69 \times 10^{-9} \frac{\alpha_0 n_n}{\mu_A^{1/2}} ,$$
 (127)

$$\mu_A = \frac{A_i A_n}{A_i + A_n} \,, \tag{128}$$

where  $\alpha_0$  is the polarizability and  $\nu_{in}$  is in s<sup>-1</sup>. Several temperature dependent ion-neutral collision frequencies from it Bailey and Balan [1996] are used in place of (126); specifically

$$\nu_{\text{O}^{+}\text{O}} = 4.45 \times 10^{-11} n(\text{O}) T^{1/2} (1.04 - 0.067 \log_{10} T)^{2},$$
(129)

where  $T = (T_{O^+} + T_n)/2$ ;

$$\nu_{\text{H}^+\text{O}} = 6.61 \times 10^{-11} n(\text{O}) T^{1/2} (1 - 0.047 \log_{10} T)^2,$$
(130)

where  $T = T_{H^+}$ ;

$$\nu_{\rm N_2^+N_2} = 5.14 \times 10^{-11} n({\rm N_2}) T^{1/2} (1.04 - 0.069 {\rm log_{10}} T)^2, \eqno(131)$$

where  $T = (T_{N_2^+} + T_n)/2$ ; and

$$\nu_{\mathrm{O_2^+O_2}} = 2.59 \times 10^{-11} n(\mathrm{O_2}) T^{1/2} (1.04 - 0.073 \log_{10} T)^2, \tag{132}$$

where  $T = (T_{O_2^+} + T_n)/2$ .

The general form of the electron-neutral heating term is

$$Q_{en} = \sum_{q} \frac{2m_e m_q}{(m_e + m_q)^2} \nu_{eq} (T_q - T_e) , \qquad (133)$$

which we rewrite as

$$Q_{en} = \sum_{q} \bar{\nu}_{eq} \left( T_q - T_e \right) , \qquad (134)$$

and  $\bar{\nu}_{eq}$  denotes the heating rate. The following rates are used in SAMI2 where the units are eV cm<sup>-3</sup> s<sup>-1</sup>. The elastic electron-neutral rates are given by [Banks and Kockarts, 1973]

$$\bar{\nu}_{eN_2} = 8.00 \times 10^{-20} n(\mathrm{N_2}) (1 - 1.2 \times 10^{-4} T_e) T_e \ , \ (135)$$

$$\bar{\nu}_{eO_2} = 5.27 \times 10^{-19} n(O_2) (1 + 3.6 \times 10^{-2} T_e^{1/2}) T_e^{1/2} ,$$
(136)

$$\bar{\nu}_{eO} = 4.8 \times 10^{-18} n(\text{O}) T_e^{1/2} ,$$
 (137)

$$\bar{\nu}_{eH} = 4.2 \times 10^{-16} n(\mathrm{H}) (1 - 1.35 \times 10^{-4} T_e) T_e^{1/2} \ . \ (138)$$

The rotational excitation rates are given by [Banks and Kockarts, 1973]

$$\bar{\nu}_{\rm N_2}^{\rm rot} = 1.33 \times 10^{-14} n({\rm N_2}) T_e^{-1/2} ,$$
 (139)

$$\bar{\nu}_{\rm O_2}^{\rm rot} = 4.67 \times 10^{-14} n({\rm O_2}) T_e^{-1/2} \ .$$
 (140)

The vibrational excitation rate is given by  $[Millward\ et\ al.,\ 1996]$ 

$$\bar{\nu}_{N_2}^{\text{vib}} = 4.33 \times 10^{-22} n(N_2) (T_n - 310)^2 \exp[0.0023 (T_e - T_n)]$$
 (141)

Finally, the electron-ion collision frequency is given by [Kelley, 1989]

$$\nu_{ei} = \left[34 + 4.18 \ln \left(T_e^3/n_e\right)\right] n_n T_e^{-3/2}$$
 (142)

where n is in cm<sup>-3</sup>,  $T_e$  is in kelvins, and  $\nu_{ei}$  is in s<sup>-1</sup>.

**Acknowledgments.** We thank P. Bernhardt, G. Bailey, and P. Richards for several helpful discussions during the development of SAMI2. This research has been supported by the Office of Naval Research under the Accelerated Research Initiative 'Ionospheric Specification and Forecasting.'

Janet G. Luhmann thanks Graham Bailey and Dwight T. Decker for their assistance in evaluating this paper.

