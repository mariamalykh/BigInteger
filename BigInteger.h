#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

const long long base = 1000000000;

class BigInteger {
 public:

  BigInteger() {}

  std::vector<long long> data() {
    return digits_;
  }

  std::vector<long long> data() const {
    return digits_;
  }


  friend bool operator==(const BigInteger& lbi, const BigInteger& rbi);
  friend bool operator<(const BigInteger& lbi, const BigInteger& rbi);
  friend bool operator<=(const BigInteger& lbi, const BigInteger& rbi);
  friend BigInteger operator*(const BigInteger& lbi, const BigInteger& rbi);
  friend BigInteger operator-(const BigInteger& lbi, const BigInteger& rbi);
  BigInteger operator-() const;
  void RemoveZeros();
  void SumPositives(const BigInteger& number);
  void DifferPositives(const BigInteger& number);
  void shift();
  friend BigInteger GCD(const BigInteger& a, const BigInteger& b);
  friend int main();

  bool isPositive() {
    return is_positive_;
  }

  bool isPositive() const {
    return is_positive_;
  }

  inline void swap(BigInteger& a, BigInteger& b) noexcept {
    std::swap(a.is_positive_, b.is_positive_);
    std::swap(a.digits_, b.digits_);
  }

  bool IsLessWithoutSign(const BigInteger& lbi, const BigInteger& rbi) {
    if (lbi.data().size() != rbi.data().size()) {
      return (lbi.data().size() < rbi.data().size());
    }

    for (int i = lbi.data().size() - 1; i >= 0; i--) {
      if (lbi.data()[i] != rbi.data()[i]) {
        return (lbi.data()[i] < rbi.data()[i]);
      }
    }
    return false;
  }

  BigInteger(long long value) {
    long long num;
    if (value < 0) {
      is_positive_ = false;
      num = value * (-1);
    } else {
      is_positive_ = true;
      num = value;
    }
    digits_.push_back(num % base);
    num /= base;
    while (num != 0) {
      digits_.push_back(num % base);
      num /= base;
    }
  }

  BigInteger(std::string str) {
    if (str.length() == 0) {
      digits_.push_back(0);
      is_positive_ = 1;
    } else {
      is_positive_ = 1;
      if (str[0] == '-') {
        is_positive_ = 0;
        str = str.substr(1);
      }
      while (str.length() % 9 != 0) {
        str = '0' + str;
      }
      for (size_t i = 0; i < str.length(); i += 9) {
        long long v = 0;
        for (size_t j = i; j < i + 9; ++j) {
          v = v * 10 + (str[j] - '0');
        }
        digits_.push_back(v);
      }
    }
    std::reverse(digits_.begin(), digits_.end());
  }

  std::string addition(long long& a) const {
    std::string res = std::to_string(a);
    while (res.length() != 9) {
      res = '0' + res;
    }
    return res;
  }

  bool isZero() const {
    return (digits_.size() == 1 && digits_[0] == 0);
  }

  void fixSign() { *this *= (2 * is_positive_ - 1); }

  std::string toString() const {
    BigInteger num(*this);
    std::string res;
    if (num.digits_.empty() || (num.digits_.size() == 1 && num.digits_[0] == 0)) {
      return "0";
    }
    if (num.digits_.size() == 1) {
      if (!num.is_positive_) {
        res += '-';
      }
      res += std::to_string(num.digits_[0]);
      return res;
    }
    bool check = 1;
    if (!num.is_positive_) {
      check = 0;
    }
    for (size_t i = 0; i < num.digits_.size() - 1; ++i) {
      long long digit = num.digits_[i];
      res = addition(digit) + res;
    }
    res = std::to_string(num.digits_[num.digits_.size() - 1]) + res;
    if (!check) {
      res = '-' + res;
    }
    return res;
  }

  BigInteger& operator+=(const BigInteger& number) {
    if (is_positive_ == number.is_positive_) {
      SumPositives(number);
      RemoveZeros();
      if (isZero()) {
        is_positive_ = 1;
      }
      return *this;
    } else {
      if (!is_positive_) {
        is_positive_ = 1;
        if (number < *this) {
          DifferPositives(number);
          is_positive_ = 0;
          RemoveZeros();
          if (isZero()) {
            is_positive_ = 1;
          }
          return *this;
        } else {
          auto copy = number;
          copy.DifferPositives(*this);
          *this = copy;
          RemoveZeros();
          if (isZero()) {
            is_positive_ = 1;
          }
          return *this;
        }
      } else {
        auto copy = number;
        copy.is_positive_ = 1;
        if (copy < *this) {
          DifferPositives(copy);
          RemoveZeros();
          if (isZero()) {
            is_positive_ = 1;
          }
          return *this;
        } else {
          copy.DifferPositives(*this);
          *this = copy;
          is_positive_ = 0;
          RemoveZeros();
          if (isZero()) {
            is_positive_ = 1;
          }
          return *this;
        }
      }
    }
  }

  BigInteger& operator-=(const BigInteger& number) {
    BigInteger copy = number;
    copy.is_positive_ = !copy.is_positive_;
    *this += copy;
    return *this;
  }

  BigInteger& operator*=(const BigInteger& number) {
    long long t = 0;
    BigInteger result = 0;
    result.digits_.resize(digits_.size() + number.digits_.size(), 0);
    for (size_t i = 0; i < digits_.size(); ++i) {
      for (size_t j = 0; j < number.digits_.size(); ++j) {
        result.digits_[i + j] += digits_[i] * number.digits_[j];
        result.digits_[i + j + 1] += result.digits_[i + j] / base;
        result.digits_[i + j] %= base;
      }
    }
    for (size_t i = 0; i < result.digits_.size(); ++i) {
      result.digits_[i] += t;
      t = result.digits_[i] / base;
      result.digits_[i] %= base;
    }
    if (t) {
      digits_.push_back(t);
    }
    if (is_positive_ != number.is_positive_) {
      result.is_positive_ = 0;
    }
    *this = result;
    RemoveZeros();

    return *this;
  }

  BigInteger& operator/=(const BigInteger& number);

  BigInteger& operator%=(const BigInteger& number) {


    BigInteger particular = *this;
    particular /= number;
    particular *= number;
    particular.RemoveZeros();
    *this -= particular;
    if (!is_positive_) {
      is_positive_ = 0;
    }
    return *this;
  }

  BigInteger& operator++() {
    *this += BigInteger(1);
    return *this;
  }

  BigInteger operator++(int) {
    BigInteger copy = *this;
    *this += BigInteger(1);
    return copy;
  }

  BigInteger& operator--() {
    *this -= BigInteger(1);
    return *this;
  }

  BigInteger operator--(int) {
    BigInteger copy = *this;
    *this -= BigInteger(1);
    return copy;
  }

  explicit operator bool() {
    return !(*this == 0);
  }

 private:
  std::vector<long long> digits_;
  bool is_positive_;
};

BigInteger& BigInteger::operator/=(const BigInteger& number) {
  BigInteger b = number;
  b.is_positive_ = true;
  BigInteger result = 0;
  BigInteger actual = 0;
  result.digits_.resize(digits_.size(), 0);
  for (long long i = static_cast<long long>(digits_.size()) - 1; i >= 0; --i) {
    actual.shift();
    actual.digits_[0] = digits_[i];
    actual.RemoveZeros();
    long long x = 0, left = 0, right = base;
    while (left <= right) {
      long long m = (left + right) / 2;
      BigInteger t = b * m;
      if (t <= actual) {
        x = m;
        left = m + 1;
      }
      else {
        right = m - 1;
      }
    }

    result.digits_[i] = x;
    actual -= b * x;
  }
  result.is_positive_ = is_positive_ == number.is_positive_;
  result.RemoveZeros();
  *this = result;

  return *this;
}

BigInteger operator+(const BigInteger& lbi, const BigInteger& rbi) {


  auto copy = lbi;
  copy += rbi;
  return copy;
}

BigInteger operator/(const BigInteger& lbi, const BigInteger& rbi) {

  auto copy = lbi;
  copy /= rbi;
  return copy;
}

BigInteger operator-(const BigInteger& lbi, const BigInteger& rbi) {


  auto copy = lbi;
  copy -= rbi;
  return copy;
}

BigInteger operator*(const BigInteger& lbi, const BigInteger& rbi) {


  auto copy = lbi;
  copy *= rbi;
  return copy;

}

BigInteger operator%(const BigInteger& lbi, const BigInteger& rbi) {


  auto copy = lbi;
  copy %= rbi;
  return copy;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& num) {
  auto copy = num;
  copy.RemoveZeros();
  std::string str = copy.toString();
  out << str;
  return out;
}

std::istream& operator>>(std::istream& in, BigInteger& num) {
  std::string number;
  in >> number;
  num = BigInteger(number);
  return in;
}

bool operator==(const BigInteger& lbi, const BigInteger& rbi) {


  if (rbi.is_positive_ != lbi.is_positive_) {
    return false;
  }
  if (rbi.data() == lbi.data()) {
    return true;
  }
  return false;
}

bool operator!=(const BigInteger& lbi, const BigInteger& rbi) {
  return !(lbi == rbi);
}

bool operator<(const BigInteger& lbi, const BigInteger& rbi) {


  if (lbi.is_positive_ == 0 && rbi.is_positive_ == 1) {
    return true;
  } else if (lbi.is_positive_ == 1 && rbi.is_positive_ == 0) {
    return false;
  } else if (lbi.is_positive_ == 1 && rbi.is_positive_ == 1) {
    if (lbi.digits_.size() != rbi.digits_.size()) return (lbi.digits_.size() < rbi.digits_.size());
    for (int i = lbi.digits_.size() - 1; i >= 0; i--) {
      if (lbi.digits_[i] != rbi.digits_[i]) return (lbi.digits_[i] < rbi.digits_[i]);
    }
  }
  return IsLessWithoutSign(rbi, lbi);
}

bool operator>(const BigInteger& lbi, const BigInteger& rbi) {
  return (lbi != rbi && !(lbi < rbi));
}

bool operator<=(const BigInteger& lbi, const BigInteger& rbi) {
  return !(lbi > rbi);
}

bool operator>=(const BigInteger& lbi, const BigInteger rbi) {
  return !(lbi < rbi);
}

BigInteger operator "" _bi(const char* string, size_t) {
  BigInteger result(string);
  return result;
}

BigInteger operator "" _bi(unsigned long long int number) {
  BigInteger result(number);
  return result;
}

BigInteger BigInteger::operator-() const {
  auto copy = *this;
  if (copy == 0) {
    return copy;
  }

  copy.is_positive_ = !copy.is_positive_;
  return copy;
}

void BigInteger::RemoveZeros() {
  while (digits_.size() > 1 && digits_[digits_.size() - 1] == 0) {
    digits_.pop_back();
  }
}

void BigInteger::SumPositives(const BigInteger& number) {
  if (is_positive_ == number.is_positive_) {
    int t = 0, s;
    int n = digits_.size();
    int m = number.digits_.size();
    if (m > n) {
      for (int i = 0; i < (m - n); ++i) {
        digits_.push_back(0);
      }
    }
    n = digits_.size();
    for (int i = 0; i < n; ++i) {
      if (i < m) {
        s = (digits_[i] + number.digits_[i]) + t;
      } else {
        s = digits_[i] + t;
      }
      t = s / base;
      digits_[i] = s % base;
    }
    if (t) {
      digits_.push_back(t);
    }
    return;
  }
}

void BigInteger::DifferPositives(const BigInteger& number) {
  int len1 = digits_.size();
  int len2 = number.digits_.size();
  int carry = 0, s;
  for (int i = 0; i < len1; ++i) {
    if (i >= len2) {
      s = digits_[i] + carry;
    } else {
      s = digits_[i] - number.digits_[i] + carry;
    }
    if (s >= 0) {
      carry = 0;
    } else {
      s += base;
      carry = -1;
    }
    digits_[i] = s;
  }
  RemoveZeros();
  return;
}

void BigInteger::shift() {
  if (this->digits_.size() == 0) {
    this->digits_.push_back(0);
    return;
  }
  this->digits_.push_back(this->digits_[this->digits_.size() - 1]);
  for (size_t i = this->digits_.size() - 2; i > 0; --i) this->digits_[i] = this->digits_[i - 1];
  this->digits_[0] = 0;
}


class Rational {
 public:

  BigInteger num() const {
    return numerator;
  }

  BigInteger den() const {
    return denominator;
  }

  friend bool operator<(const Rational& lhs, const Rational& rhs);
  friend bool operator==(const Rational& lhs, const Rational& rhs);
  Rational operator-() const;
  friend BigInteger GCD(BigInteger& a, BigInteger& b);
  friend int main();

  Rational(const BigInteger& a) {
    denominator = 1;
    is_positive_ = (a >= 0);
    numerator = a * (2 * is_positive_ - 1);
  }
  Rational(const long long& a) {
    denominator = 1;
    is_positive_ = (a >= 0);
    numerator = a * (2 * is_positive_ - 1);
  }

  Rational() {}

  bool isPositive() const {
    return is_positive_;
  }

  Rational& operator+=(const Rational& number) {
    if (!is_positive_) {
      numerator = -numerator * number.denominator + denominator * number.numerator * (2 * number.is_positive_ - 1);
    } else {
      numerator = numerator * number.denominator + denominator * number.numerator * (2 * number.is_positive_ - 1);
    }
    denominator *= number.denominator;
    is_positive_ = (numerator.isPositive() == denominator.isPositive());

    fixSign();

    Reduction(numerator, denominator);

    if (numerator.isZero()) {
      is_positive_ = true;
    }
    return *this;
  }

  Rational& operator*=(const Rational& number) {
    numerator *= number.numerator;
    denominator *= number.denominator;
    is_positive_ = (isPositive() == number.isPositive());

    fixSign();

    Reduction(numerator, denominator);

    if (numerator.isZero()) {
      is_positive_ = true;
    }
    return *this;
  }

  Rational& operator/=(const Rational& number) {
    numerator *= number.denominator;
    denominator *= number.numerator;
    is_positive_ = (isPositive() == number.isPositive());

    fixSign();

    Reduction(numerator, denominator);

    if (numerator.isZero()) {
      is_positive_ = true;
    }
    return *this;
  }

  Rational& operator-=(const Rational& number) {
    *this += (-(number));
    return *this;
  }

  Rational& operator++() {
    *this += Rational(1);
    return *this;
  }

  Rational operator++(int) {
    Rational copy = *this;
    *this += Rational(1);
    return copy;
  }

  Rational& operator--() {
    *this += Rational(-1);
    return *this;
  }

  Rational operator--(int) {
    Rational copy = *this;
    *this -= Rational(1);
    return copy;
  }

  void fixSign() {
    if (numerator < 0) {
      numerator *= -1;
    }

    if (denominator < 0) {
      denominator *= -1;
    }
  }

  std::string toString() const {
    std::string result;

    if (numerator.isZero()) {
      return numerator.toString();
    }

    if (!is_positive_) {
      result += '-';
    }

    result += numerator.toString();
    if (denominator != 1) {
      result += "/" + denominator.toString();
    }
    return result;
  }

  explicit operator double() {
    return std::stod(asDecimal(30));
  }

  std::string asDecimal(size_t precision = 0) const {
    BigInteger fraction = numerator;
    BigInteger num = numerator;
    BigInteger den = denominator;
    std::string mantiss;
    std::string antie;
    if (!is_positive_) {
      mantiss = "-";
    }
    mantiss += (fraction / den).toString();
    fraction %= den;
    while (antie.size() < precision) {
      fraction *= 10;
      antie += (fraction / den).toString();
      fraction %= den;
    }
    return mantiss + "." + antie;
  }

  void Reduction(BigInteger& numerator, BigInteger& denominator) {
    bool flag = 1;
    if (numerator.isPositive() != denominator.isPositive()) {
      flag = 0;
    }
    if (numerator.isPositive() == 0) {
      numerator = -numerator;
    }
    if (denominator.isPositive() == 0) {
      denominator = -denominator;
    }
    BigInteger gcd = GCD(numerator, denominator);
    numerator /= gcd;
    denominator /= gcd;
    if (flag == 0 && numerator != 0) {
      numerator = -numerator;
    }
  }
 private:
  BigInteger numerator;
  BigInteger denominator;
  bool is_positive_;
};

bool operator<(const Rational& lhs, const Rational& rhs) {
  if (lhs.is_positive_ != rhs.is_positive_) {
    return lhs.is_positive_ < rhs.is_positive_;
  }

  BigInteger copy_lhs = lhs.numerator * rhs.denominator;
  BigInteger copy_rhs = rhs.numerator * lhs.denominator;

  if (lhs.is_positive_) {
    return copy_lhs < copy_rhs;
  }
  return  copy_lhs > copy_rhs;
}

bool operator>(const Rational& lhs, const Rational& rhs) {
  return rhs < lhs;
}

bool operator==(const Rational& lhs, const Rational& rhs) {
  if (lhs.is_positive_ != rhs.is_positive_) {
    return false;
  }

  BigInteger lhs_copy = lhs.numerator * rhs.denominator;
  lhs_copy.fixSign();
  BigInteger rhs_copy = rhs.numerator * lhs.denominator;
  rhs_copy.fixSign();

  std::cerr << "result: " << (rhs_copy == lhs_copy) << '\n';

  return rhs_copy == lhs_copy;
}

bool operator!=(const Rational& lhs, const Rational& rhs) {
  return !(lhs == rhs);
}

bool operator<=(const Rational& lhs, const Rational& rhs) {
  return !(lhs > rhs);
}

bool operator>=(const Rational& lhs, const Rational& rhs) {
  return !(lhs < rhs);
}

Rational Rational::operator-() const {
  if (numerator == 0) {
    return *this;
  }

  auto copy = *this;
  copy.is_positive_ = !copy.is_positive_;
  return copy;
}

BigInteger GCD(const BigInteger& a, const BigInteger& b) {
  if (a % b == 0) {
    return b;
  }
  if (b % a == 0) {
    return a;
  }
  if (a > b) {
    return GCD(a % b, b);
  }
  return GCD(a, b % a);
}


Rational operator+(const Rational& lbi, const Rational& rbi) {
  auto copy = lbi;
  copy += rbi;
  return copy;
}

Rational operator-(const Rational& lbi, const Rational& rbi) {
  auto copy = lbi;
  copy -= rbi;
  return copy;
}

Rational operator/(const Rational& lbi, const Rational& rbi) {
  auto copy = lbi;
  copy /= rbi;
  return copy;
}

Rational operator*(const Rational& lbi, const Rational& rbi) {
  auto copy = lbi;
  copy *= rbi;
  return copy;
}

int main() {
  BigInteger a = BigInteger(1);
  BigInteger b = BigInteger(2);
  std::cout << b - a;
}


