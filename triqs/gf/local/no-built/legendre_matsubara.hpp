
class GF_Bloc_ImFreq;
class GF_Bloc_ImTime;
class GF_Bloc_ImLegendre;

void legendre_matsubara_direct (GF_Bloc_ImLegendre const & Gl, GF_Bloc_ImFreq & Gw);
void legendre_matsubara_inverse (GF_Bloc_ImFreq const & Gw, GF_Bloc_ImLegendre & Gl);

void legendre_matsubara_direct (GF_Bloc_ImLegendre const & Gl, GF_Bloc_ImTime & Gt);
void legendre_matsubara_inverse (GF_Bloc_ImTime const & Gt, GF_Bloc_ImLegendre & Gl);

