/// SortCard.cpp : This file contains the 'main' function. Program execution begins and ends there.
// RJM 6/2/24

#include "cardlib.h"

#include <string>
#include <iostream>
#include <cstdlib>
#include <ctime> // for time()
using namespace std;


const int maxCard = 20;
int compareCount = 0;
int swapCount = 0;
int maxRecursionLevel = 0;

aCard thePack[maxCard];

void printPack(std::string mess)
{
    std::cout << mess << " : ";
    for (int ct = 0; ct < maxCard; ct++)
    {
        std::cout << cardToStr(thePack[ct]) << " ";
        if (ct != maxCard - 1) std::cout << ", ";

    }
}

int compareCards(aCard c1, aCard c2){
    compareCount++;
    if (c1.cardSuit < c2.cardSuit) {
        return -1;
    } else if (c1.cardSuit > c2.cardSuit) {
        return 1;
    } else { // Suits are equal, compare by value
        if (c1.cardVal < c2.cardVal) {
            return -1;
        } else if (c1.cardVal > c2.cardVal) {
            return 1;
        } else {
            return 0; // Cards are equal
        }
    }
}

void swapCards(aCard& c1, aCard& c2){
    swapCount++;
    aCard temp = c1;
    c1 = c2;
    c2 = temp;
}

void bubbleSort(aCard a[], int n){
    compareCount = 0;
    swapCount = 0;
    for (int i = 0; i < n-1; i++){
        for (int j = 0; j < n-i-1; j++){
            if (compareCards(a[j], a[j+1]) == 1){
                swapCards(a[j], a[j+1]);
            }
        }
    }
}

int partition(aCard a[], int low, int high, int& compareCount, int& swapCount) {
    aCard pivot = a[high];
    int i = low - 1;

    for (int j = low; j <= high - 1; j++) {
        compareCount++;
        if (compareCards(a[j], pivot) == -1) {
            i++;
            swapCards(a[i], a[j]);
            swapCount++;
        }
    }

    swapCards(a[i + 1], a[high]);
    swapCount++;
    return i + 1;
}


void quickSort(aCard a[], int low, int high, int& compareCount, int& swapCount, int& maxRecursionLevel) {
    if (low < high) {
        int partitionIndex = partition(a, low, high, compareCount, swapCount);

        if (partitionIndex - low > maxRecursionLevel) {
            maxRecursionLevel = partitionIndex - low;
        }

        quickSort(a, low, partitionIndex - 1, compareCount, swapCount, maxRecursionLevel);
        quickSort(a, partitionIndex + 1, high, compareCount, swapCount, maxRecursionLevel);
    }
}



int main() {
    cout << "Card Sorting!\n";

    for (int ct = 0; ct < maxCard; ct++)
        thePack[ct] = getCard("32009512");  // generate a random pack
    printPack("\nUnsorted\n");



    bubbleSort(thePack, maxCard); // Pass thePack and maxCard as arguments to bubbleSort

    //quickSort(thePack, 0, maxCard - 1, compareCount, swapCount, maxRecursionLevel);

    printPack("\n\nSorted");
    cout << "\n\nComparisons: " << compareCount << " Swaps: " << swapCount << " Max Recursion Level: " << maxRecursionLevel << endl;

    return 0;
}

