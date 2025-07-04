#pragma once


template<T>
int index_of_mode_value(std::vector<T> array)
{
  int number = array[0];
  int mode = number;
  int count = 0;
  int index = 0;
  int countMode = 1;

  for (int i=1; i<array.size(); i++)
  {
    if (array[i] == number) 
    { // count occurrences of the current number
      ++count;
    }
    else
    { // now this is a different number
      if (count > countMode) 
      {
        countMode = count; // mode is the biggest occurance
        mode = number;
        index = i;
      }
      count = 1; // reset count for the new number
      number = array[i];
    }
  }

  return index;
}

int index_of_mode_value(std::vector<int> array)
{
  int count = array[0];
  int index = 0;

  for(int i=1; i<array.size(); i++) {
    if( array[i] > count ){
      count = array[i];
      index = 1;
    }
  }




