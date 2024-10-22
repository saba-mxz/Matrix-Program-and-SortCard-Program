// RJM's file for getting a random card
// 6/2/24

#include "cardlib.h"
using namespace std;

string cardToStr(aCard c) {
    string res = "  ";

    if (c.cardVal == 1) res[0] = 'A';
    else if (c.cardVal < 10) res[0] = '0' + c.cardVal; // convert int to string and take first char (0..9
    else if (c.cardVal == 10) res[0] = 'T';
    else if (c.cardVal == 11) res[0] = 'J';
    else if (c.cardVal == 12) res[0] = 'Q';
    else if (c.cardVal == 13) res[0] = 'K';
    else res[0] = 'X';

    if (c.cardSuit == hearts) res[1] = 'H';
    else if (c.cardSuit == clubs) res[1] = 'C';
    else if (c.cardSuit == diamonds) res[1] = 'D';
    else if(c.cardSuit == spades) res[1] = 'S';

    return res;

}


aCard getCard(string stdno)
{  // function returns a card  - using student number stdno (8 numerical chars)
    aCard ans;
    string rcardstr = stdno.substr(rand() % 8, 1) + stdno.substr(rand() % 8, 1);    // two random characters from stdno
    int rcard = stoi(rcardstr) % 52;  // get integer in range 0..51
    string res = "  ";
    ans.cardVal = 1 + rcard % 13; // give cardVal 1..13
    ans.cardSuit = static_cast<Suit>(rcard / 13); // and for suit
    return ans;
}

//
// Created by Saba Mirza on 22/10/2024.
//
