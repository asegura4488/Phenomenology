#include <iostream>
#include <fstream>
#include <string>

using namespace std; 

class particle
{
  public:
    int id;
    int status;
    int mother1, mother2;
    int color1, color2;
    double px, py, pz, p0;
    double m;
    double t;
    double spin;
    double pt;
    double eta;
};

int main(int argc, char *argv[])
{
  string card, s1;

  int n_particles, i1;
  double d1, d2, d3, d4;

  particle all_particles[1000];

  if (argc != 3) cout << "Please specify an input file and an output file! Quitting..." << endl;

  ifstream in_file(argv[1], ios::in);
  ofstream out_file(argv[2], ios::out);

  do
  {
    getline(in_file, card); 
    out_file << card << endl;   

    s1 = card.substr(0,7);
  }
  while (s1 != "</init>");

  while (in_file >> card)
  {
    if (card == "<event>")
    {
      out_file << card << endl;
      in_file >> n_particles >> i1 >> d1 >> d2 >> d3 >> d4;
      out_file << n_particles << " " << i1 << " " << d1 << " " << d2 << " " << d3 << " " << d4 << endl; 
  
      for (int i = 0; i < n_particles; i++)
      {
        in_file >> all_particles[i].id     >> all_particles[i].status >> all_particles[i].mother1 >> all_particles[i].mother2
                >> all_particles[i].color1 >> all_particles[i].color2 >> all_particles[i].px      >> all_particles[i].py 
                >> all_particles[i].pz     >> all_particles[i].p0     >> all_particles[i].m       >> all_particles[i].t        >> all_particles[i].spin;     

      }

      for (int i = 0; i < n_particles; i++)
      {
        if (all_particles[i].id == 21 && all_particles[i].status > 0 && all_particles[all_particles[i].mother1 - 1].id != 6000001 && all_particles[all_particles[i].mother2 - 1].id != 6000001)
        {
          all_particles[i].color2 = all_particles[i].color1;
        }
      }

      for (int i = 0; i < n_particles; i++)
      {
        out_file <<  all_particles[i].id << " " << all_particles[i].status << " " << all_particles[i].mother1 << " " << all_particles[i].mother2 << " " 
                << all_particles[i].color1 << " " << all_particles[i].color2 << " " << all_particles[i].px      << " " << all_particles[i].py << " " 
                << all_particles[i].pz     << " " << all_particles[i].p0     << " " << all_particles[i].m       << " " << all_particles[i].t        << " " << all_particles[i].spin << endl;
      }
    }
    else if (card == "</LesHouchesEvents>")
    {
      out_file << card << endl;
      break;
    }
    else if (card == "</event>")
    {
      out_file << card << endl;
    }
    else
    {
      cout << card << endl;
    }
  }

  in_file.close();
  out_file.close();

  return 0;
}