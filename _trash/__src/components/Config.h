#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include "Window.h"


using std::vector;
using std::string;

//typedef unordered_map<string, double>;

class Config
{
private:
    vector<Window> windows;
    
public:
    Config();
    Config(vector<Window> given);
    Config(vector<Window> given, Window without);
    void add(Window win);
    double totalWeight();
    //double adjustedWeight();
    const vector<Window> getVector() const;
    bool overlapsWindow(Window win);
    int size() const;
    vector<Window>::iterator begin();
    vector<Window>::iterator end();
    friend std::ostream& operator<<(std::ostream& os, const Config& C);
    void print();
    Window getItem(int i);
    double Z();
    void sort();
    string toString();

    bool operator== (Config &I2);
    bool operator!= (Config &I2);   
    friend bool operator== (const Config& wA, const Config& wB);

    friend Config sorted(const Config & c);

    static Config fromString(string s);

};

namespace std {
template <>
  struct hash<Config>
  {
    std::size_t operator()(const Config& k) const
    {
      using std::size_t;
      using std::hash;
      using std::string;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:
      string s = "";
      for(auto w : k.getVector())
        s += w.toString();
      cout << "here22 = " << k << endl;
      cout << (hash<string>()(s)) << endl;
      return (hash<string>()(s));
    }
  };
}


#endif