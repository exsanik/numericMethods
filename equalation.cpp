class equalation
{
private:
	int equalPower;
	double *equalat;
	double borders[2];
	double eps = 0.0000001;

	double _max(const double a, const double b) { (a > b) ? a : b; }
	double f1(double x, double h = 0.00001) { return (f(x + h) - f(x - h)) / (2 * h); } //Первая производная
	double f2(double x, double h = 0.00001) { return (f(x + h) - 2 * f(x) + f(x - h)) / (h * h); } //Вторая производная

public:
	equalation()
	{
		equalPower = 0;
		equalat = NULL;
		borders[2] = { 0 };
	}
	equalation(int power, const double * eq)
	{
		equalPower = power + 1;
		equalat = new double[equalPower];
		for (int i = 0; i < equalPower; i++)
			equalat[i] = eq[i];
		bordersEq(*this);
	}
	equalation(equalation& obj)
	{
		equalPower = obj.equalPower;
		borders[0] = obj.borders[0];
		borders[1] = obj.borders[1];
		for (int i = 0; i < equalPower; i++)
			equalat[i] = obj.equalat[i];
	}
	~equalation()
	{
		delete[] equalat;
	}
	double* bordersEq(equalation& obj)
	{
		double max = abs(obj.equalat[1]);
		for (int i = 2; i < obj.equalPower; i++)
			if (abs(obj.equalat[i]) > max)
				max = abs(obj.equalat[i]);
		obj.borders[1] = 1 + max / abs(obj.equalat[0]);
		obj.borders[0] = -obj.borders[1];
		return obj.borders;
	}
	void setBorders(const double a, const double b)
	{
		borders[0] = a;
		borders[1] = b;
	}
	void setEqualation(int power, const double * eq)
	{
		equalPower = power + 1;
		equalat = new double[equalPower];
		for (int i = 0; i < equalPower; i++)
			equalat[i] = eq[i];
		bordersEq(*this);
	}
	double f(double x)
	{
		double f = 0;
		for (int i = 0; i < equalPower; i++)
			f += equalat[i] * pow(x, equalPower - i - 1);
		return f;
	}

	double iteration()
	{
		double m = f1(borders[0]), M = f1(borders[1]);
		if (m > M) swap(m, M);
		double lamb = 1 / M, q = 1 - (m / M);
		double x = 0, y;
		int iter = 0;
		do
		{
			y = x;
			x = x - f(x)*lamb;
			iter++;
			cout << "Step = " << iter << " x = " << myDouble << x << endl;
		} while (abs(x - y) > abs((1 - q) * eps / q));
		return x;
	}

	double tangents() 
	{
		double x1 = borders[1] - f(borders[1]) / f1(borders[1]);
		double x0 = borders[1];
		int iter = 0;
		while (abs(x0 - x1) > eps) {
			x0 = x1;
			x1 = x1 - f(x1) / f1(x1);
			iter++;
			cout << "Step = " << iter << " x = " << myDouble << x1 << endl;
		}
		return x1;
	}
	double *hornerCoefficient(double coef)
	{
		double *res = new double[equalPower];
		res[0] = equalat[0];
		for (int i = 1; i < equalPower; i++)
		{
			double ch = coef * res[i - 1] + equalat[i];
			res[i] = (abs(ch) < eps)? 0 : ch;
		}
		equalPower--;
		delete[] equalat;
		equalat = res;
		return equalat;
	}

	void solveTangents()
	{
		while(equalPower > 1)
		{
			cout << "Tangents method " << endl;
			double res1 = this->tangents();
			cout << "Result tangents method = " << res1 << endl;
			this->hornerCoefficient(res1);
			this->bordersEq(*this);
		}
	}
	void solveIter()
	{
		while (equalPower > 1)
		{
			cout << "Iteration method " << endl;
			double res1 = this->iteration();
			cout << "Result iteration method = " << res1 << endl;
			this->hornerCoefficient(res1);
			this->bordersEq(*this);
		}
	}
};
