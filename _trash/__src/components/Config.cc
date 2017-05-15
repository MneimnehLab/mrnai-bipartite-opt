#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include "Window.h"
#include "Config.h"


using std::vector;
using std::string;


Config::Config()
{
}

Config::Config(vector<Window> given)
{
    for(auto win : given)
        windows.push_back(win);
}

Config::Config(vector<Window> given, Window without)
{
    for(auto win : given)
    {
        if(win != without)
            windows.push_back(win);
    }
}

void Config::add(Window win)
{

    windows.push_back(win);
}


double Config::totalWeight()
{
    double weight = 0;
    for(auto win : windows)
        weight += win.weight;

    return weight;
}


const vector<Window> Config::getVector() const
{
    return windows;
}


bool Config::overlapsWindow(Window win)
{
    bool overlaps = false;
    for(auto our : windows)
        if(our.overlaps(win))
            overlaps = true;

    return overlaps;
}

int Config::size() const
{
    return windows.size();
}


vector<Window>::iterator Config::begin()
{
    return windows.begin();
}

vector<Window>::iterator Config::end()
{
    return windows.end();
}

void Config::print()
{

    cout << "[" ;
    for(int i=0;i < (int)windows.size();i++)
    {
        cout << windows[i];
        if(i < (int)windows.size() - 1) 
            cout << ", " ;
    }
    cout << "]" << endl;
}

Window Config::getItem(int i)
{
    return windows[i];
}


std::ostream& operator<<(std::ostream& os, const Config& C)
{
    int i=0;
    os << "[" ;
    for(;i< (int)C.windows.size();i++)
    {
        os << C.windows.at(i);
        if(i < (int)C.windows.size() - 1) 
            os << ", ";
    }
    os << "]";
    return os;
}


/*
double Config::Z()
{
    return exp((double)totalWeight());
}
*/



bool Config::operator== (Config &I2)
{
    if(windows.size() != I2.windows.size())
        return false;
    else
    {
        for(int i=0;i<windows.size();i++)
        {
            if(windows[i] != I2.windows[i])
                return false;
        }
    }
    return true;
}

bool operator== (const Config& wA, const Config& I2)
{
    auto windows = wA.windows;
    if(windows.size() != I2.windows.size())
        return false;
    else
    {
        for(int i=0;i<windows.size();i++)
        {
            if(windows[i] != I2.windows[i])
                return false;
        }
    }
    return true;
}

bool Config::operator!= (Config &I2)
{
    return !(*this == I2);
}

bool myfunction2 (Window i, Window j) { return (i<j); }


void Config::sort()
{
    std::sort(windows.begin(), windows.end(), myfunction2);
}

Config  sorted(const Config & c)
{
    vector<Window> v(c.windows);
    std::sort(v.begin(), v.end(), myfunction2);
    return Config(v);
}

string Config::toString()
{
    stringstream ss;
    ss.str("");
    for(auto w : windows)
        ss << w.toString();

    return ss.str();
      
}

Config Config::fromString(string line)
{
    Config c;

    for(int i=0;i<line.length();i++)
        if ( !isdigit(line.at(i)) )
            line[i] = ' ';

    std::stringstream ss;
    ss.str(line);
    int l,i,j,w1,w2;
    //int count = 0;
    while(ss >> l >> i >> j >> w1 >> w2)
    {
        Window Window(l, i, j, w1, w2);

        if(w1 == w2)
            c.add(Window);
    }

    return c;
}


