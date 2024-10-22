#pragma once

// header file for cardlib
// RJM 06/01/24


#include <string>

enum Suit { hearts, clubs, diamonds, spades };      // define suits

struct aCard {                          // defines a card
    int cardVal;                        // number 1..13
    Suit cardSuit;                      // suit
};

std::string cardToStr(aCard c);			// declare function to represent a card as a two character string

aCard getCard(std::string stdno);        // declares function to get a card from stdno string

int compareCards(aCard c1, aCard c2);