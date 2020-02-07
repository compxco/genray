#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cctype>
#include <vector>
#include <algorithm>

using namespace std;

/////////////////
// files stuff //
/////////////////

void file_copy(string from, string to)
{
  remove(to.c_str());
  rename(from.c_str(),to.c_str());
}

//////////////////
// string stuff //
//////////////////

// from pos to first space, skip initial spaces
string str_extract(string& s, int pos)
{
  int l = s.length();
  string t = "";
  int i;
  for( i=pos; i<l&&isspace(s[i]); i++)
    ;
  for( ; i<l&&!isspace(s[i]); i++ )
    t += s[i];
  return t;
}

// identifier from pos
string str_extract_name(string& s, int pos)
{
  int l = s.length();
  string t = "";
  int i;
  for( i=pos; i<l&&isspace(s[i]); i++)
    ;
  for( ; i<l&&(isalnum(s[i])||s[i]=='_'); i++ )
    t += s[i];
  return t;
}

// from pos to end-1
string str_extract(string& s, int pos, int end)
{
  int l = s.length();
  string t = "";
  for( int i=pos; i<end; i++ )
    t += s[i];
  return t;
}

// from pos to end of line, skip initial spaces
string str_extract2(string& s, int pos)
{
  int i;
  int l = s.length();
  for( i=pos; i<l&&isspace(s[i]); i++)
    ;
  string t="";
  for( ; i<l; i++ )
    t += s[i];
  return t;
}

string str_to_lower(string s)
{
  string r="";
  int l=s.length();
  for( int i=0; i<l; i++ )
    r += tolower(s[i]);
  return r;
}

string space_line(int l)
{
  string r="";
  for( int i=0; i<l; i++ )
    r += ' ';
  return r;
}

string int_to_str(int n)
{
  char buf[100];
  sprintf(buf,"%d",n);
  return string(buf);
}

/////////////////
// lines stuff //
/////////////////

/* [cC]<some_text> */
bool is_comment(string& s)
{
  return (s[0]=='c' || s[0]=='C');
}

/* [cC]<type> <some_text> */
bool is_comment(string& s, string& type)
{
  if( !is_comment(s) )
    return false;
  string t = str_extract(s,1);
  return t==type;
}

/* <some_text>!**<type> */
bool is_uncomment(string& s, string& type, int& pos)
{
  pos = s.find("!**");
  string t;
  if( pos<0 )
    return false;
  t = str_extract(s,pos+3);
  return t==type;
}

/* check empty lines */
bool is_empty_my(string& s)
{
  int l = s.length();
  for( int i=0; i<l; i++ )
    if( !isspace(s[i]) )
      return false;
  return true;
}

/* check line is continuation */
bool is_cont(string& s)
{
  if( is_comment(s) || is_empty_my(s) || s.length()<6 )
    return false;
  return str_extract(s,0,5)=="     " && s[5]!=' ';
}

// makes right alignment in uncommented line
string align_line(string& s)
{
  // remove all tabulations
  string t="";
  int l = s.length();
  for( int i=0; i<l; i++ )
    if( s[i]!='\t' )
      t += s[i];
    else
      t += "      ";
  // special cases
  if( is_comment(t) || is_empty_my(t) )
    return t;
  if( is_cont(t) )
    {
      string r = "     ";
      r += "+"+str_extract2(t,6);
      return r;
    }
  // align line
  string prefix = str_extract(t,0,6);
  string rest = str_extract2(t,6);
  return prefix+rest;
}

string remove_end_comment(string s)
{
  string t = "";
  int l=s.length();
  bool state = true;
  for( int i=0; i<l; i++ )
    {
      if( s[i]=='!' && state )
	return t;
      if( s[i]=='"' )
	state = !state;
      t += s[i];
    }
  return t;
}

///////////////////
// actions stuff //
///////////////////

void do_comments_of_type(string type, string fname)
{
  ifstream f(fname.c_str());
  ofstream g("tmp");
  
  string line;
  char buf[256];
  
  f.getline(buf,255);
  int pos;

  while( f )
    {
      line = buf;
      if( !is_uncomment(line,type,pos) )
	g << line << endl;
      else
	{
	  line = str_extract(line,0,pos);
	  g << "C" << type << " " << line << endl;
	}
      f.getline(buf,255);
    }
 
  g.close();
  f.close();
  file_copy("tmp",fname);
}

void undo_comments_of_type(string type, string fname)
{
  ifstream f(fname.c_str());
  ofstream g("tmp");
  
  string line;
  char buf[256];
  
  f.getline(buf,255);
  int pos;

  while( f )
    {
      line = buf;
      if( !is_comment(line,type) )
	g << line << endl;
      else
	{
	  line = str_extract(line,2+type.length(),line.length());
	  g << line << "!**" << type << endl;
	}
      f.getline(buf,255);
    }
 
  g.close();
  f.close();
  file_copy("tmp",fname);
}

void remove_comments_of_type(string type, string fname)
{
  ifstream f(fname.c_str());
  ofstream g("tmp");
  
  string line;
  char buf[256];
  
  f.getline(buf,255);

  while( f )
    {
      line = buf;
      if( !is_comment(line,type) )
	g << line << endl;
      f.getline(buf,255);
    }
 
  g.close();
  f.close();
  file_copy("tmp",fname);
}

void remove_comment_marks_of_type(string type, string fname)
{
  ifstream f(fname.c_str());
  ofstream g("tmp");
  
  string line;
  char buf[256];
  
  f.getline(buf,255);
  int pos;

  while( f )
    {
      line = buf;
      if( !is_comment(line,type) )
	g << line << endl;
      else
	{
	  line = str_extract(line,2+type.length(),line.length());
	  g << line << endl;
	}
      f.getline(buf,255);
    }
 
  g.close();
  f.close();
  file_copy("tmp",fname);
}

void remove_all_comments(string fname)
{
  ifstream f(fname.c_str());
  ofstream g("tmp");
  
  string line;
  char buf[256];
  
  f.getline(buf,255);

  while( f )
    {
      line = buf;
      if( !is_comment(line) )
	g << line << endl;
      f.getline(buf,255);
    }

  g.close();
  f.close();
  file_copy("tmp",fname);
}

void clear_buf(vector<string>& big_line, ostream& g, string pat, string type)
{
  string prefix = "C";
  prefix += type+" ";
  bool c = false;
  int l = big_line.size();
  if( l==0 )
    return;
  for( int i=0; i<l; i++ )
    {
      int pos = big_line[i].find(pat);
      c = c || (pos>=0);
    }
  string p;
  if( c )
    p = prefix;
  else
    p = "";
  for( int i=0; i<l; i++ )
    g << p << big_line[i] << endl;
  big_line.clear();
}

void do_comments_of_pattern(string pat, string type, string fname)
{
  ifstream f(fname.c_str());
  ofstream g("tmp");
  
  string line;
  char buf[256];
  vector<string> big_line;
  
  f.getline(buf,255);

  while( f )
    {
      line = buf;
      if( !is_comment(line) )
	{
	  if( !is_cont(line) )
	    clear_buf(big_line,g,pat,type);
	  big_line.push_back(line);
	}
      else
	{
	  clear_buf(big_line,g,pat,type);
	  g << line << endl;
	}
      f.getline(buf,255);
    }
  clear_buf(big_line,g,pat,type);

  g.close();
  f.close();
  file_copy("tmp",fname);
}

void align_all_lines(string fname)
{
  ifstream f(fname.c_str());
  ofstream g("tmp");
  
  string line;
  char buf[256];
  
  f.getline(buf,255);

  while( f )
    {
      line = buf;
      g << align_line(line) << endl;
      f.getline(buf,255);
    }

  g.close();
  f.close();
  file_copy("tmp",fname);
}

enum { _program, _subroutine, _function };

////////////////
// main stuff //
////////////////

void print_info(int key=-1)
{
  cerr << "usage:" << endl;
  cerr << " cmm -R sources              : remove all commented lines  " << endl;
  cerr << " cmm -r type sources         : remove all lines commented by type" << endl;
  cerr << " cmm -c type sources         : make comments of type" << endl;
  cerr << " cmm -u type sources         : clear all comments of type" << endl;
  cerr << " cmm -e type sources         : remove comment marks of type" << endl;
  cerr << " cmm -p type pattern sources : make comments for all lines with pattern" << endl;
  cerr << " cmm -a sources              : align all uncommented lines" << endl;
  if( key!=-1 )
    exit(key);
}

int main(int argc, char** argv)
{
  if( argc==1 )
    print_info(1);

  string action = argv[1];
  if( action=="-R" )
    {
      if( argc<3 )
	print_info(2);
      for( int i=2; i<argc; i++ )
	remove_all_comments(argv[i]);
      return 0;
    }
  if( action=="-r" )
    {
      if( argc<4 )
	print_info(2);
      string type = argv[2];
      for( int i=3; i<argc; i++ )
	remove_comments_of_type(type,argv[i]);
      return 0;
    }
  if( action=="-c" )
    {
      if( argc<4 )
	print_info(2);
      string type = argv[2];
      for( int i=3; i<argc; i++ )
	do_comments_of_type(type,argv[i]);
      return 0;
    }
  if( action=="-u" )
    {
      if( argc<4 )
	print_info(2);
      string type = argv[2];
      for( int i=3; i<argc; i++ )
	undo_comments_of_type(type,argv[i]);
      return 0;
    }
  if( action=="-p" )
    {
      if( argc<5 )
	print_info(2);
      string type = argv[2];
      string pat = argv[3];
      for( int i=4; i<argc; i++ )
	do_comments_of_pattern(pat,type,argv[i]);
      return 0;
    }
  if( action=="-a" )
    {
      if( argc<3 )
	print_info(2);
      for( int i=2; i<argc; i++ )
	align_all_lines(argv[i]);
      return 0;
    }
  if( action=="-e" )
    {
      if( argc<4 )
	print_info(2);
      string type = argv[2];
      for( int i=3; i<argc; i++ )
	remove_comment_marks_of_type(type,argv[i]);
      return 0;
    }
  print_info();
  return 2;
}
