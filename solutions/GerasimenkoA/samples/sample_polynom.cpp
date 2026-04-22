// ННГУ, ИИТММ, Курс "Алгоритмы и структуры данных"
//
// Copyright (c) Пинежанин Е.С.
//
// Тестирование работы с полиномами

#include <iostream>

#include "tpolynom.h"

int main()
{
  TPolynom p1;
  p1.SetPolynom("3x^2y^4 + 4.5yz^2");

  TPolynom p2;
  p2.SetPolynom("3x^2y^4 + 4.5yz^2");

  TPolynom p3 = p1 + p2;

  std::cout << p3 << " = " << p3.Calculate(1.0, 2.5, 4.7) << "\n";

  return 0;
}
