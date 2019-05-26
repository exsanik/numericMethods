class Interpolation
{
private:
	function<double(double)> fun;
	double a, b, h;

	vector<pair<double, double>> table;
	vector<pair<double, double>> lagrange;
	vector<pair<double, double>> newtone;

public:
	Interpolation(function<double(double)> f, double a, double b) : fun(f), a(a), b(b)
	{
		h = (b - a) / 10;
		buildTable();
		buildLagrange();
		buildNewton();
	}

	double Lagrange(double x)
	{
		double P = 0;
		for (int j = 0; j < table.size(); j++)
		{
			double Q = 1;
			for (int i = 0; i < table.size(); i++)
			{
				if (i != j)
					Q *= (x - table[i].first) / (table[j].first - table[i].first);
			}
			P += table[j].second * Q;
		}
		return P;
	}

	double Newton2(double x)
	{
		vector<double> deltay;
		for (int j = 0; j < table.size(); j++)
			deltay.push_back(table[j].second);

		double P = table[table.size() - 1].second, S = 1;
		double t = (x - table[table.size()-1].first) / h;
		for (int i = table.size() - 1, k = 1; i > 0; i--, k++)
		{	
			for (int j = 0; j < i; j++)
			{
				deltay[j] = deltay[j + 1] - deltay[j];
				if (abs(deltay[0] - deltay[1]) < 0.0000001 && abs(deltay[1] - deltay[2]) < 0.0000001)
					break;
			}

			S *= (t + k-1) / k;
			P += S * deltay[i-1];
		}
		return P;
	}

	double Newton1(double x)
	{
		double sum = table[0].second;
		vector<double> d;
		vector<double> deltay;
		for (int j = 0; j < table.size(); j++)
			deltay.push_back(table[j].second);
		for (int i = 0; i < table.size(); i++)
		{
			for (int j = 0; j < deltay.size() - i - 1; j++)
				deltay[j] = deltay[j + 1] - deltay[j];
			d.push_back(deltay[0]);
			if (abs(deltay[0] - deltay[1]) < 0.0000001 && abs(deltay[1] - deltay[2]) < 0.0000001)
				break;
		}

		double t = (x - table[0].first) / h;
		double S = 1;
		for (int j = 0; j < d.size(); j++)
		{
			S *= (t - j) / (j + 1);
			sum += d[j] * S;
		}

		return sum;
	}

	void outputCoords()
	{
		for (int i = 0; i < table.size(); i++)
		{
			cout << "Function: x = " << table[i].first << "\ty = " << table[i].second << endl;
			cout << "Lagrange: x = " << lagrange[i].first << "\ty = " << lagrange[i].second << endl;
			cout << "Newtone: x = " << newtone[i].first << "\ty = " << newtone[i].second << endl;
			cout << endl;
		}
	}

	vector<pair<double, double>> getLagrange() { return lagrange; }
	vector<pair<double, double>> getNewtone() { return lagrange; }
	vector<pair<double, double>> getCoords() { return table; }


protected:

	void buildTable()
	{
		for (double step = a; step <= b; step += h)
			table.emplace_back(pair<double, double>(step, fun(step)));
	}

	void buildLagrange()
	{
		for (auto &coord : table)
		{
			pair<double, double> lagCoord;
			lagCoord.first = coord.first;
			lagCoord.second = Lagrange(lagCoord.first);
			lagrange.emplace_back(lagCoord);
		}
	}

	void buildNewton()
	{
		for (auto &coord : table)
		{
			pair<double, double> newtCoord;
			newtCoord.first = coord.first;
			newtCoord.second = Newton2(newtCoord.first);
			newtone.emplace_back(newtCoord);
		}
	}
};
