class Poisson2D {
private:
  static const int dim = 2;

public:
  class DiffusionCoefficient {
  public:
    double
    value(const unsigned int i, const unsigned int j) {
      return i + j;
    }
  };

  class ForcingTerm {
  public:
    double
    value(const unsigned int i, const unsigned int j) {
      return i - j;
    }
  };
};
