// ННГУ, ИИТММ, Курс "Алгоритмы и структуры данных"
//
// tcomputecentermodel.cpp
//
// Copyright (c) Пинежанин Е.С.

#include <tpolynom.h>
#include <stdexcept>  
#include <istream>
using namespace std;

TPolynom::TPolynom() {};

int packPowers(int x, int y, int z) {
	return x * 10000 + y * 100 + z;
}

int getPowerX(int powers) { return powers / 10000; }
int getPowerY(int powers) { return (powers / 100) % 100; }
int getPowerZ(int powers) { return powers % 100; }

bool TPolynom::IsCorrect() const {
	if (data.empty()) return true;

	list<int> powers;
	for (auto i = data.begin(); i != data.end(); ++i) {
		const TTerm& term = *i;

		if (term.coeff == 0.0) return false;

		int px = getPowerX(term.powers);
		int py = getPowerY(term.powers);
		int pz = getPowerZ(term.powers);
		if (px < 0 || py < 0 || pz < 0) return false;

		if (find(powers.begin(), powers.end(), term.powers) != powers.end()) return false;

		powers.push_back(term.powers);
	}
	return true;
}

void TPolynom::sort_polynom() {
	if (data.empty()) return;

	data.sort([](const TTerm& a, const TTerm& b) { return a.powers > b.powers; });

	auto iter = data.begin();
	while (iter != data.end()) {
		auto next = iter;
		next++;
		//если степени одинаковые, то пред. += след. и удаляем следующий
		if (next != data.end() && iter->powers == next->powers) { 
			iter->coeff += next->coeff;
			data.erase(next);
		} else {
			iter++;
		}
	}
	data.remove_if([](const TTerm& tmp) {return tmp.coeff == 0.0; });
}

//Было принято решение вынести парсинг в отдельную ф-ию, тк SetPolynom и Адд используют это.. (дабы не писать два раза)
static TTerm parsing(const string& monom) {
	string s;
	for (char c : monom) if (!isspace(c)) s += c;
	if (s.empty()) throw out_of_range("Empty monom");

	TTerm term;
	term.coeff = 1.0;
	size_t i = 0;
	if (i < s.length() && (isdigit(s[i]) || s[i] == '.')) {
		size_t start = i;
		while (i < s.length() && (isdigit(s[i]) || s[i] == '.')) i++;
		term.coeff = stod(s.substr(start, i - start)); //stod - string to double :)
	}
	int px = 0, py = 0, pz = 0;
	while (i < s.length() && (s[i] == 'x' || s[i] == 'y' || s[i] == 'z')) {
		char var = s[i];
		i++;
		int exp = 1;
		if (i < s.length() && s[i] == '^') {
			i++;
			if (i >= s.length() || !isdigit(s[i])) throw out_of_range("Error '^'");
			size_t startExp = i;
			while (i < s.length() && isdigit(s[i])) i++;
			exp = stoi(s.substr(startExp, i - startExp)); //stoi = string to int :)
		}
		if (var == 'x') px = exp;
		else if (var == 'y') py = exp;
		else if (var == 'z') pz = exp;
	}
	term.powers = packPowers(px, py, pz);
	return term;
}

void TPolynom::SetPolynom(const std::string& polynom) {
	data.clear();
	string s;
	for (char c : polynom) {
		if (!isspace(c)) s += c;
	}
	if (s.empty()) throw out_of_range("Empty string");

	size_t i = 0;
	bool first = true;
	while (i < s.length()) {
		double sign = 1.0;
		if (!first && (s[i] == '+' || s[i] == '-')) {
			sign = (s[i] == '+') ? 1.0 : -1.0;
			i++;
		}
		else if (first && s[i] == '-') {
			sign = -1.0;
			i++;
		}
		first = false;

		if (i >= s.length()) throw out_of_range("Unexpected end of string");

		// Находим конец монома (до следующего + или -)
		size_t start = i;
		while (i < s.length() && s[i] != '+' && s[i] != '-') i++;
		if (start == i) throw out_of_range("Empty monom");

		string monomStr = s.substr(start, i - start);
		TTerm term = parsing(monomStr);
		term.coeff *= sign;
		data.push_back(term);
	}
	sort_polynom();
}

TPolynom TPolynom::operator+(const TPolynom &polynom) const {
	//Просто добавляем новые члены в конец и "сортируем"

	TPolynom result = *this;
	result.data.insert(result.data.end(), polynom.data.begin(), polynom.data.end());
	result.sort_polynom();
	return result;
}

TPolynom TPolynom::operator-(const TPolynom &polynom) const {
	//Добавляем новые члены с противоположным знаком и "сортируем"

	TPolynom result = *this;
	for (const TTerm& term : polynom.data) {
		TTerm tmp = term;
		tmp.coeff = -tmp.coeff;
		result.data.push_back(tmp);
	}
	result.sort_polynom();
	return result;
}

TPolynom TPolynom::operator*(const TPolynom &polynom) const {
	TPolynom result;
	for (const TTerm& a : data) {
		for (const TTerm& b : polynom.data) {
			TTerm term;
			term.coeff = a.coeff * b.coeff;
			int px = getPowerX(a.powers) + getPowerX(b.powers);
			int py = getPowerY(a.powers) + getPowerY(b.powers);
			int pz = getPowerZ(a.powers) + getPowerZ(b.powers);
			term.powers = packPowers(px, py, pz);
			result.data.push_back(term);
		}
	}
	result.sort_polynom();
	return result;
}

TPolynom TPolynom::operator*(double coeff) const {
	if (coeff == 0.0) return TPolynom();

	TPolynom result = *this;
	for (TTerm& term : result.data) {
		term.coeff *= coeff;
	}
	return result;
}

void TPolynom::Add(const std::string& monom) {
	if (monom.empty()) throw out_of_range("Empty monom");
	TTerm term = parsing(monom);
	data.push_back(term);
	sort_polynom();
}

void TPolynom::Delete(size_t pos) {
	if (pos >= data.size()) throw out_of_range("Invalid position");
	auto iter = data.begin();
	for (size_t i = 0; i < pos; ++i) ++iter;
	data.erase(iter);
}

double TPolynom::Calculate(double x, double y, double z) const {
	double result = 0.0;
	for (const TTerm& term : data) {
		int px = getPowerX(term.powers);
		int py = getPowerY(term.powers);
		int pz = getPowerZ(term.powers);
		result += term.coeff * pow(x, px) * pow(y, py) * pow(z, pz);
	}
	return result;
}

std::istream& operator>>(std::istream& istr, TPolynom& polynom) {
	string line;
	istr >> ws; // убираем пробелы/переводы строк перед чтением (с <istream>)
	getline(istr, line);
	polynom.SetPolynom(line);
	return istr;
}

std::ostream& operator<<(std::ostream& ostr, const TPolynom& polynom) {
	if (polynom.data.empty()) {
		ostr << "0";
		return ostr;
	}

	bool first = true;
	for (const TTerm& term : polynom.data) {
		if (!first && term.coeff > 0) ostr << "+";
		first = false;

		int px = getPowerX(term.powers);
		int py = getPowerY(term.powers);
		int pz = getPowerZ(term.powers);

		if (term.coeff != 1 && term.coeff != -1) {
			ostr << term.coeff;
		}
		else if (term.coeff == -1) {
			ostr << "-";
		}

		if (px > 0) ostr << "x^" << px;
		if (py > 0) ostr << "y^" << py;
		if (pz > 0) ostr << "z^" << pz;

		if (px == 0 && py == 0 && pz == 0) {
			if (term.coeff == 1) ostr << "1";
			else if (term.coeff == -1) ostr << "-1";
		}
	}
	return ostr;
}
