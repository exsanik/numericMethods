class Integral
{
private:
	double a, b, h;
	int n = 5; 
	double eps = 0.0001;
	function<double(double)> f;
public:
	Integral(double a, double b, function<double(double)> f) : a(a), b(b), f(f) { }
	double trapeziumMethod()
	{
		double part = n, step = (b - a) / part;
		double result0 = 0, result1 = 0;
		do
		{
			result1 = result0;
			result0 = (f(a) + f(a + (part*step))) / 2;
			for (double i = (a + step); i < (a + (part*step)); i += step)
				result0 += f(i);
			result0 *= step;
			
			part *= 2;
			step = (b - a) / part;
		} while ((1.0/3*abs(result0 - result1)) > eps);
		return result1;
	}

	double simpsonMethod()
	{
		double part = n, step = (b - a) / (part);
		double result0 = 0, result1 = 0;
		do
		{
			result1 = result0;
			result0 = f(a) + f(a + (part*step));
			double even = 0;
			for (double i = (a + step*2); i < (a + (part*step)); i += 2*step)
				even += f(i);
			even *= 2;
			result0 += even;
			double odd = 0;
			for (double i = (a + step); i < (a + (part*step)); i += 2*step)
				odd += f(i);
			odd *= 4;
			result0 += odd;
			result0 *= (step/3);

			part *= 2;
			step = (b - a) / part;
		} while (1.0/15*abs(result0 - result1) > eps);

		return result1;
	}
};
