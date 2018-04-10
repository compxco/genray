#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
  ifstream f[2];
  f[0].open(argv[1]);
  f[1].open(argv[2]);

  int cf = 0; // current file
  
  string buf;
  getline(f[0],buf);
  while(f[0])
    {
      if( buf.substr(0,18) != "CMPIINSERTPOSITION" )
	cout << buf << endl;
      else
	cf = 1-cf;
      getline(f[cf],buf);
    }

  f[0].close();
  f[1].close();
}

