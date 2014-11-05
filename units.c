
extern NFWParZ Halo;
#ifndef POTENTIAL_HEADER_INCLUDED
/*constants in alternative GADGET units
 * Mass: 10^10Msun *
 * Length: kpc *
 * vel:   km/s *  */
#define G 43007.1
#define HUBBLE0 0.073  //h=0.73, H0=h*0.1 (km/s/kpc)

/* virial definitions */
#define VIR_TH 0
#define VIR_C200 1
#define VIR_B200 2

void set_units(float MassInMsunh, float LengthInKpch, float VelInKms)
{
  float G0=43007.1;
  float H0=0.1;
  Units.MassInMsunh=MassInMsunh;
  Units.LengthInKpch=LengthInKpch;
  Units.VelInKms=VelInKms;
  Units.Const.G=G0*VelInKms*VelInKms*LengthInKpch/MassInMsunh;
  Units.Const.H0=H0*LengthInKpch/MassInMsunh;
}