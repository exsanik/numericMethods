typedef vector<double> eq;
typedef vector<eq> slau;

class Slau
{
private:
	double eps = 0.00001;
	slau a;
	eq b;
	eq x;
public:
	Slau(slau coef, eq res) : a(coef), b(res) { x.resize(a.size()); }
	Slau(const Slau &obj) : a(obj.a), b(obj.b), x(obj.x) {}
	eq solveIteration(Slau s)
	{
		cout << "Iteration method: " << endl;
		moveMaxToMainDiag(s);
		bool exit = true;
		int step = 1;
		while (exit)
		{
			cout << "Step: " << step++ << ":\t";
			for (int i = 0; i < s.a.size(); i++)
			{
				double x1 = s.x[i];
				s.x[i] = fn(i);
				cout << "x" << i + 1 << " = " << fixed << s.x[i] << '\t';
				if (abs(s.x[i] - x1) < eps)
				{
					exit = false;
					break;
				}
			}
			cout << endl;
		}
		outputRes(s);
		return s.x;
	}

	eq gaus(Slau s)
	{
		cout << "Gaus method: " << endl;
		for (int i = 0; i < s.a.size() - 1; i++)
		{
			for (int k = 0; k < s.a.size(); k++)
			{
				double div = s.a[k][i];
				for (int j = i; j < s.a[i].size(); j++)
					s.a[k][j] /= div;
				s.b[k] /= div;
			}

			for (int j = i+1; j < s.a[i].size(); j++)
				subStr(s, i, j);
			cout << "Step: " << i + 1 << endl;
			outputMatr(s);
			cout << endl;
		}
		reverseRoot(s);
		outputRes(s);
		return s.x;
	}

	eq mainElements(Slau s)
	{
		cout << "Main elements method: " << endl;
		cout << "Step: " << 0 << endl;
		outputMatr(s);
		cout << endl;
		for (int i = 0; i < s.a.size() - 1; i++)
		{
			int mInd = maxElement(s, i);
			eq mcoef(s.a.size(), 1);
			double main = -s.a[mInd][i];
			for (int j = i; j < s.a.size(); j++)
			{
				if (j == mInd)
					continue;
				mcoef[j] = s.a[j][i] / main;
			}
		
			for (int j = i; j < s.a.size(); j++)
			{
				if (j == mInd)
					continue;
				for (int k = i; k < s.a[j].size(); k++)
					s.a[j][k] += mcoef[j] * s.a[mInd][k];
				s.b[j] += mcoef[j] * s.b[mInd];
			}
			cout << "Step: " << i + 1 << endl;
			outputMatr(s);
			cout << endl;
		}
		reverseRoot(s);

		outputRes(s);
		return s.x;
	}
protected:
	double fn(int xN)
	{
		x[xN] = 0;
		double res = b[xN];
		for (int i = 0; i < a[xN].size(); i++)
			res -= x[i] * a[xN][i];
		res /= a[xN][xN];
		x[xN] = res;
		return x[xN];
	}

	void moveMaxToMainDiag(Slau &s)
	{
		for (int i = 0; i < s.a.size(); i++)
		{
			auto m = max_element(s.a[i].begin(), s.a[i].end());
			auto it = s.a[i].begin() + i;
			if (m != it)
				swap(*it, *m);
		}
	}

	void swapStr(Slau &s, int col)
	{
		int ind = 0;
		for (int i = 1; i < s.a.size(); i++)
			if (abs(s.a[i][col]) > abs(s.a[ind][col]))
				ind = i;
		if (ind != 0)
			swap(s.a[0], s.a[ind]);
	}

	void subStr(Slau &s, int s1, int s2)
	{
		for (int i = 0; i < s.a[s1].size(); i++)
			s.a[s2][i] -= s.a[s1][i];
		s.b[s2] -= s.b[s1];
	}

	double maxElement(Slau &s, int col)
	{
		double ind = 0;
		for (int i = 1; i < s.a.size(); i++)
			if (abs(s.a[i][col]) > abs(s.a[ind][col]))
				ind = i;
		return ind;
	}

	void reverseRoot(Slau &s)
	{
		s.x[s.x.size() - 1] = s.b[s.b.size() - 1] / s.a[s.a.size() - 1][s.a[0].size() - 1];
		for (int i = s.x.size() - 2; i >= 0; i--)
		{
			double sum = 0;
			for (int j = s.x.size() - 1; j > i; j--)
				sum += s.a[i][j] * s.x[j];

			s.x[i] = (s.b[i] - sum) / s.a[i][i];
		}
	}

	void outputMatr(Slau &s)
	{
		int i = 0;
		for (auto &el : s.a)
		{
			for (auto &e : el)
				cout << e << '\t';
			cout << s.b[i++] << endl;
		}
	}

	void outputRes(Slau &s)
	{
		cout << "Result: ";
		for (int i = 0; i < s.x.size(); i++)
			cout << "x" << i + 1 << " = " << s.x[i] << '\t';
		cout << endl<< endl;
	}
};
