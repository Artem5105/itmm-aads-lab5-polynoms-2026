// ННГУ, ИИТММ, Курс "Алгоритмы и структуры данных"
//
// polynom.h
//
// Copyright (c) Пинежанин Е.С.

#ifndef __TPolynom_H__
#define __TPolynom_H__

#include <list>
#include <string>

struct TTerm
{
  double coeff; 
  int powers; // вида ABCDEF(***), где A - степень x, B - степень y и C - степень z. 
};

class TPolynom {
private:
  std::list<TTerm> data;
  TTerm parsing(const std::string& monom);
  std::list<TTerm> ParsePolynomString(const std::string& s);

public:
  TPolynom();

  bool IsCorrect() const; // возвращает коррекность полинома (false - если на вход подали некорректное выражение)

  void CollectPolynom(); // используется для оптимизации полинома (коэфф. одинаковых степеней складываются (вычитаются))

  void SetPolynom(const std::string& polynom);

  //(У операторов убрал &, так вроде правильнее, если не +=, -= и тд)
  TPolynom operator+(const TPolynom &polynom) const; 
  TPolynom operator-(const TPolynom &polynom) const;

  TPolynom operator*(const TPolynom &polynom) const;
  TPolynom operator*(double coeff) const;
  friend TPolynom operator*(double coeff, const TPolynom& p);

  bool operator==(const TPolynom& other) const;
  bool operator!=(const TPolynom& other) const;

  void Add(const std::string& monom); // Добавление монома (убрал const - добавление монома изменяет полином, а конст не позволит)
  void Delete(size_t pos); // Удаление монома

  double Calculate(double x, double y, double z) const;

  friend std::istream& operator>>(std::istream& istr, TPolynom& polynom);
  friend std::ostream& operator<<(std::ostream& ostr, const TPolynom& polynom);
};

#endif
