#include <vector>
#include <string>
#include <sstream>
#include "Window.h"

Window::Window(int l1_, int l2_, int i_, int j_, int w1_, int w2_) 
        : l1(l1_), l2(l2_), i(i_), j(j_), w1(w1_), w2(w2_), weight(0), intrEnergy(0) {
            // cout << "hello" << endl;
        }

Window::Window(int l1_, int l2_, int i_, int j_, int w1_, int w2_, double weight_, double intrEnergy_) 
        : l1(l1_), l2(l2_), i(i_), j(j_), w1(w1_), w2(w2_), weight(weight_), intrEnergy(intrEnergy_) {
            // cout << "bye " << l1 << ", " << l2 << endl;
        }

Window::Window(int l1, int l2, int i_, int j_, int w1_, int w2_, double weight_) 
        : l1(l1), l2(l2), i(i_), j(j_), w1(w1_), w2(w2_), weight(weight_), 
          intrEnergy(0) 
{
}

Window::Window(std::string line) 
{
    //format: l i j w1 w2 delta intrEnergy
    /*
    std::stringstream ss;
    ss.str(line);
    
    ss >> l >> i >> j >> w1 >> w2;

    if(ss.rdbuf()->in_avail() == 0)
        weight = 0;
    else
        ss >> weight;
    */
    throw std::runtime_error("Constructor Window::Window(std::string line) is no longer handled in paired regions.");

}

// int Window::numRNAs = 0;

bool Window::lessThan(const Window& win1, const Window& win)
{
    if(win1.l1 < win.l1)
        return true;
    if(win1.l1 > win.l1)
        return false;
    
    if(win1.l2 < win.l2)
        return true;
    if(win1.l2 > win.l2)
        return false;
    

    if(win1.i < win.i)
        return true;
    if(win1.i > win.i)
        return false;
    
    if(win1.j < win.j)
        return true;
    if(win1.j > win.j)
        return false;
    
    if(win1.w1 < win.w1)
        return true;
    if(win1.w1 > win.w1)
        return false;
    
    if(win1.w2 < win.w2)
        return true;
    
    return false;
}

bool Window::operator< (Window &win)
{
    return Window::lessThan(*this, win);
}


bool Window::operator== (Window &I2)
{
    if(l1 == I2.l1 && l2 == I2.l2 && i == I2.i && j == I2.j && w1 == I2.w1 && w2 == I2.w2)
        return true;
    
    return false;
}

bool operator== (const Window& wA, const Window& wB)
{
    if(wB.l1 == wA.l1 && wB.l2 == wA.l2 && wB.i == wA.i && wB.j == wA.j && wB.w1 == wA.w1 && wB.w2 == wA.w2)
        return true;
    
    return false;

}

bool operator!= (const Window& wA, const Window& wB)
{
    return !(wA == wB);
}

bool operator< (const Window& wA, const Window& wB)
{
    // return wA < wB;
    return Window::lessThan(wA, wB);;
}

bool Window::operator!= (Window &I2)
{
    return !(*this == I2);
}

std::ostream& operator<<(std::ostream& os, const Window& I)
{
    os << "(" << I.l1 << ", " << I.l2 << ", " << I.i << ", " << I.j << ", " << I.w1 << ", " << I.w2 << ")";
    return os;
}

string Window::toString() const
{
    stringstream ss;
    ss.str("");
    ss << "(" << l1 << ", " << l2 << ", " << i << ", " << j << ", " << w1 << ", " << w2 << ")";
    return ss.str();
}

string Window::toStringW() const
{
    stringstream ss;
    ss.str("");
    ss << "(" << l1 << ", " << l2 << ", " << i << ", " << j << ", " << w1 << ", " << w2 << ") : " << weight;
    return ss.str();
}


std::vector<Window>* Window::createWindowSetFromPyString(std::string line)
{
    std::vector<Window> * ints = new std::vector<Window>();

    for(int i=0;i<line.length();i++)
        if ( !isdigit(line.at(i)) )
            line[i] = ' ';

    std::stringstream ss;
    ss.str(line);
    int l1,l2,i,j,w1,w2;
    //int count = 0;
    while(ss >> l1 >> l2 >> i >> j >> w1 >> w2)
    {
        Window Window(l1, l2, i, j, w1, w2);

        if(w1 == w2)
            ints->push_back(Window);
    }

    return ints;
}

bool Window::overlaps(Window win)
{
    // level here refers to the level between RNAs, not the RNAs themselves
    // A is this window, B is param window
    int AEvenLeft = i - w1;
    int AEvenRight = i;
    int AOddLeft = j - w2;
    int AOddRight = j;

    int BEvenLeft = win.i - win.w1;
    int BEvenRight = win.i;
    int BOddLeft = win.j - win.w2;
    int BOddRight = win.j;

    int AOdd = l2;
    int AEven = l1;

    int BOdd = win.l2;
    int BEven = win.l1;

    if(AOdd != BOdd and AEven != BEven)
        return false;


    // if both on same level
    if(AEven == BEven and AOdd == BOdd) // equiv. both evens same and both odds same
    {
        
        // We will consider it an over lap if both windows are immediately adjacent, i.e., 
        if((AEvenLeft == BEvenRight+1 and AOddLeft == BOddRight+1) or (AEvenRight == BEvenLeft-1 and AOddRight == BOddLeft-1))
            return true;

        return !( (AEvenLeft > BEvenRight and AOddLeft > BOddRight) or (AEvenRight < BEvenLeft and AOddRight < BOddLeft) );
    }
    // this window is above param window
    else if( AEven == BEven )
        return !(AEvenLeft > BEvenRight or AEvenRight < BEvenLeft) ;
    
    // this window is below param window
    else if( AOdd == BOdd )
        return !(AOddLeft > BOddRight or AOddRight < BOddLeft)  ;      
        

    

    return false;

}

 
void Window::prettyPrint()
{
    printf("Interaction b/w RNA %d & %d : [%d,%d] & [%d,%d]  -- weight = %f\n",
            l1, l2, i-w1, i, j-w2, j, weight );
}



// Window::Window(int l_, int i_, int j_, int w1_, int w2_) 
//         : l(l_), i(i_), j(j_), w1(w1_), w2(w2_), weight(0), intrEnergy(0) {}

// Window::Window(int l_, int i_, int j_, int w1_, int w2_, double weight_, double intrEnergy_) 
//         : l(l_), i(i_), j(j_), w1(w1_), w2(w2_), weight(weight_), intrEnergy(intrEnergy_) 
// {
//     rna1 = l;
//     rna2 = l+1;
// }

// Window::Window(int rna1, int rna2, int i_, int j_, int w1_, int w2_, double weight_) 
//         : rna1(rna1), rna2(rna2), i(i_), j(j_), w1(w1_), w2(w2_), weight(weight_), 
//           intrEnergy(0) 
// {
//     l = -1;
// }

// Window::Window(std::string line) 
// {
//     //format: l i j w1 w2 delta intrEnergy
//     /*
//     std::stringstream ss;
//     ss.str(line);
    
//     ss >> l >> i >> j >> w1 >> w2;

//     if(ss.rdbuf()->in_avail() == 0)
//         weight = 0;
//     else
//         ss >> weight;
//     */
//     throw std::runtime_error("Constructor Window::Window(std::string line) is no longer handled in paired regions.");

// }

// bool Window::operator< (Window &win)
// {
//     //Order matters!!!
//     if(l < win.l)
//         return true;
//     if(l > win.l)
//         return false;
    
//     if(i < win.i)
//         return true;
//     if(i > win.i)
//         return false;
    
//     if(j < win.j)
//         return true;
//     if(j > win.j)
//         return false;
    
//     if(w1 < win.w1)
//         return true;
//     if(w1 > win.w1)
//         return false;
    
//     if(w2 < win.w2)
//         return true;
    
//     return false;

// }

// bool Window::operator== (Window &I2)
// {
//     if(l == I2.l && i == I2.i && j == I2.j && w1 == I2.w1 && w2 == I2.w2)
//         return true;
    
//     return false;
// }

// bool operator== (const Window& wA, const Window& wB)
// {
//     if(wB.l == wA.l && wB.i == wA.i && wB.j == wA.j && wB.w1 == wA.w1 && wB.w2 == wA.w2)
//         return true;
    
//     return false;

// }

// bool operator!= (const Window& wA, const Window& wB)
// {
//     return !(wA == wB);
// }

// bool operator< (const Window& wA, const Window& wB)
// {
//     //Order matters!!!
//     auto win = wB;
//     if(wA.l < win.l)
//         return true;
//     if(wA.l > win.l)
//         return false;
    
//     if(wA.i < win.i)
//         return true;
//     if(wA.i > win.i)
//         return false;
    
//     if(wA.j < win.j)
//         return true;
//     if(wA.j > win.j)
//         return false;
    
//     if(wA.w1 < win.w1)
//         return true;
//     if(wA.w1 > win.w1)
//         return false;
    
//     if(wA.w2 < win.w2)
//         return true;
    
//     return false;
// }

// bool Window::operator!= (Window &I2)
// {
//     return !(*this == I2);
// }

// std::ostream& operator<<(std::ostream& os, const Window& I)
// {
//     // os << "(" << I.l << ", " << I.i << ", " << I.j << ", " << I.w1 << ", " << I.w2 << ")";
//     os << "(" << I.rna1 << ", " << I.rna2  << ", " << I.i << ", " << I.j << ", " << I.w1 << ", " << I.w2 << ")";
//     return os;
// }

// string Window::toString() const
// {
//     stringstream ss;
//     ss.str("");
//     ss << "(" << l << ", " << i << ", " << j << ", " << w1 << ", " << w2 << ")";
//     return ss.str();
// }

// string Window::toStringW() const
// {
//     stringstream ss;
//     ss.str("");
//     ss << "(" << l << ", " << i << ", " << j << ", " << w1 << ", " << w2 << ") : " << weight;
//     return ss.str();
// }


// std::vector<Window>* Window::createWindowSetFromPyString(std::string line)
// {
//     std::vector<Window> * ints = new std::vector<Window>();

//     for(int i=0;i<line.length();i++)
//         if ( !isdigit(line.at(i)) )
//             line[i] = ' ';

//     std::stringstream ss;
//     ss.str(line);
//     int l,i,j,w1,w2;
//     //int count = 0;
//     while(ss >> l >> i >> j >> w1 >> w2)
//     {
//         Window Window(l, i, j, w1, w2);

//         if(w1 == w2)
//             ints->push_back(Window);
//     }

//     return ints;
// }

// bool Window::overlaps(Window win)
// {
//     int ATopLeft = i - w1;
//     int ATopRight = i;
//     int ABotLeft = j - w2;
//     int ABotRight = j;

//     int BTopLeft = win.i - win.w1;
//     int BTopRight = win.i;
//     int BBotLeft = win.j - win.w2;
//     int BBotRight = win.j;

//     if(l == win.l)
//     {
//         //return !( (ATopLeft > BTopRight and ABotLeft > BBotRight) or (ATopRight < BTopLeft and ABotRight < BBotLeft) );

//         // We will consider it an over lap if both windows are immediately adjacent, i.e., 
//         // (ATopLeft == BTopRight+1 and ABotLeft == BBotRight+1) or (ATopRight == BTopLeft-1 and ABotRight == BBothLeft-1)
//         if((ATopLeft == BTopRight+1 and ABotLeft == BBotRight+1) or (ATopRight == BTopLeft-1 and ABotRight == BBotLeft-1))
//             return true;

//         return !( (ATopLeft > BTopRight and ABotLeft > BBotRight) or (ATopRight < BTopLeft and ABotRight < BBotLeft) );
//     }
    
//     else if(l == win.l - 1)
//         return !(ABotLeft > BTopRight or ABotRight < BTopLeft)  ;      
    
//     else if(l == win.l + 1)
//         return !(ATopLeft > BBotRight or ATopRight < BBotLeft) ;

//     return false;

// }
