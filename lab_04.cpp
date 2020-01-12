#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

class SimplexMethod {
	// строка таблицы
	struct Row {
		vector<double> a; // коэффициенты при x
		double b; // правая часть
	};

	int n; // число переменных
	int m; // число ограничений
	vector<Row> table; // симплекс таблица
	vector<double> c; // коэффициенты оптимизируемой функции
	vector<int> variables; // все переменные
	vector<int> basis; // базисные переменные
	vector<double> deltas; // дельты
	vector<double> solvation; // решение

	void CalculateDeltas(); // вычисление дельт

	int GetArgMinDelta(); // вычисление номера минимальной дельты
	int GetArgMaxDelta(); // вычисление номера максимальной дельты

	void InitialVariables(int resistanseNumber = -1); // инициализация переменных и базиса

	void MakeNewBasis(double pivot, int index, int jindex); // задание нового базисного элемента
	int MaxNegativeB(); // максимальная по модулю отрицательная b

	int CheckForIntegerSolvation(); // проверка на то, что решение целое
	int GetF(vector<int> solve); // получение значения целевой функции по вектору
	double GetF(vector<double> solve); // получение значения целевой функции по вектору
	void PrintVector(vector<int> vec); // вывод вектора
	void PrintVector(vector<double> vec); // вывод вектора
	int GetMaxFloatIndex(); // получение индекса максимального дробного числа

	vector<double> GetRestriction(int index, double &res);

	bool isSolve; // есть ли решение
	void MakeBasis(); // создание базиса
	int isBasis(int i); // проверка на наличие базиса в столбце
public:
	SimplexMethod(); // конструктор по умолчанию со всеми данными

	void SetRestriction(vector<double> restriction, double res); // добавление условия

	void Read(); // ввод значений
	void Print(); // вывод таблицы

	void Solve(int max); // решение ОЗЛП
	void RemoveNegativeB(); // удаление отрицательных b

	void IntegerSolve(); // поиск целочисленного решения
};

// конструктор по умолчанию со всеми данными
SimplexMethod::SimplexMethod() {
	// задаём значения для конкрутной задачи
	n = 3;
	m = 3; // увеличиваем кол-во условий, если нужно

	table = vector<Row>(m, { vector<double>(n), 0 });
	c = vector<double>(n);

	c = { 2, 6, 7 };

	table[0].a = { 3, 1, 1 };
	table[0].b = { 3 };

	table[1].a = { 1, 2, 0 };
	table[1].b = { 8 };

	table[2].a = { 0, 0.5, 2 };
	table[2].b = { 1 };
	solvation = vector<double>(n);

	InitialVariables(); // инициализируем переменные
	RemoveNegativeB();

	isSolve = false;
}

// добавление условия
void SimplexMethod::SetRestriction(vector<double> restriction, double res) {
	m++; // увеличиваем кол-во условий

	table.push_back({ restriction, res }); // добавляем условие в таблицу

	InitialVariables(m - 1); // инициализируем переменные
	RemoveNegativeB(); // удаляем отрицательные b

	isSolve = false; // пока не решено
}

 // проверка на наличие базиса в столбце
int SimplexMethod::isBasis(int i) {
	int zeroCount = 0; // кол-во нулей в столбце
	int basisIndex = 0; // индекс не нулевого элемента

	// идём по столбцу
	for (int j = 0; j < m; j++) {
		// если это ноль
		if (table[j].a[i] == 0)
			zeroCount++; // увеличиваем кол-во нулей
		// иначе
		else
			basisIndex = j; // запоминаем индекс
	}

	// если нули все, кроме одного, то это базис
	if (zeroCount == m - 1)
		return basisIndex;

	return -1;
}

// создание базиса
void SimplexMethod::MakeBasis() {
	// ижём по всем перемнным
	for (int i = 0; i < variables.size(); i++) {
		int basisInd = isBasis(i); // проверяем, является ли эта переменная базисной

		// если является
		if (basisInd != -1) {
			basis[basisInd] = i; // запоминаем её

			MakeNewBasis(table[basisInd].a[i], basisInd, i); // делаем на этом месте базис
		}
	}

	// идём по базисным переменным
	for (int j = 0; j < basis.size(); j++) {
		// если переменная не запонена
		if (basis[j] == -1) {
			int iIndex = 0; // индекс нового базисного элемента

			// находим ненулевой элемент в строке и запоминаем индекс
			for (int i = 0; i < variables.size(); i++) {
				if (table[j].a[i] != 0) {
					iIndex = i;
					break;
				}
			}

			MakeNewBasis(table[j].a[iIndex], j, iIndex); // делаем базис
			basis[j] = iIndex; // запоминаем
		}
	}
} 

// инициализация переменных и базиса
void SimplexMethod::InitialVariables(int resistanseNumber) {
	variables.clear(); // очищаем переменные
	basis = vector<int>(m, -1);

	c.erase(c.begin() + table[0].a.size(), c.end());
	
	// добавляем переменные
	for (int i = 0; i < n + m; i++)
		variables.push_back(i);

	int size = table[0].a.size();

	for (int i = 0; i < m; i++) {
		c.push_back(0); // добавляем нули в функцию

		// добавляем коэффициенты для переменных с коэффициентом 1, если они стоят на главной диагонали, иначе с нулём - если не передан номер условия
		// если передан номер условия, то проставляем 1 только в номере условия
		for (int j = 0; j < (m + n) - size; j++)
			table[i].a.push_back(resistanseNumber == -1 ? (i == j) : j == resistanseNumber);
	}

	MakeBasis(); // создаём базис
}


// поиск макисмальной отрицательной b
int SimplexMethod::MaxNegativeB() {
	int imax = -1;

	// ищем максимальный отрицательный элемент
	for (int i = 1; i < m; i++) {
		if (table[i].b < 0 && (imax == -1 || table[i].b < table[imax].b))
			imax = i;
	}

	return imax; // возвращаем максимум
}

// устранение отрицательной правой части
void SimplexMethod::RemoveNegativeB() {
	int imax = MaxNegativeB(); // индекс максимального по модулю отрицательного элемента

	// пока если отрицательные элементы
	while (imax != -1) {
		int jmax = 0; // индекс новой базисной переменной

					  // идём по столбцу и ищем максимальный по модул. элемент
		for (int j = 1; j < m; j++) {
			if (fabs(table[imax].a[j]) > fabs(table[imax].a[jmax]))
				jmax = j;
		}

		basis[imax] = jmax;	// запоминаем индекс новой базисной переменной
		MakeNewBasis(table[imax].a[jmax], imax, jmax); // делаем этот элемент базисным

		imax = MaxNegativeB(); // находим новый максимальный по модулю элемент в правой части
	}
}

// создание новой базисной переменной на месте index, jindex
void SimplexMethod::MakeNewBasis(double pivot, int index, int jindex) {
	// делим строку на элемент
	for (size_t i = 0; i < table[index].a.size(); i++)
		table[index].a[i] /= pivot;

	table[index].b /= pivot;

	// вычитаем из всех остальных строк эту строку, умноженную на элемент в столбце jmax
	for (int i = 0; i < m; i++) {
		if (i == index)
			continue;

		double value = table[i].a[jindex];

		for (size_t j = 0; j < table[i].a.size(); j++)
			table[i].a[j] -= table[index].a[j] * value;

		table[i].b -= table[index].b * value;
	}
}

// ввод значений
void SimplexMethod::Read() {
	cout << "Enter function coefficients (c): ";
	c = vector<double>(n); // создаём вектор коэффициентов

						   // считываем коэффициенты оптимизируемой функции
	for (int i = 0; i < n; i++)
		cin >> c[i];

	cout << "Enter restrictions coefficients:" << endl;

	// считываем коэффициенты ограничений
	for (int i = 0; i < m; i++) {
		cout << "Enter restriction " << (i + 1) << ": ";

		for (int j = 0; j < n; j++)
			cin >> table[i].a[j];

		cin >> table[i].b;
	}

	variables.clear(); // очищаем переменные

					   // добавляем переменные
	for (int i = 0; i < n; i++)
		variables.push_back(i);

	for (int i = 0; i < m; i++) {
		c.push_back(0); // добавляем нули в функцию
		variables.push_back(n + i); // добавляем доп переменные
		basis.push_back(n + i); // делаем их базисными

								// добавляем коэффициенты для переменных с коэффициентом 1, если они стоят на главной диагонали, иначе с нулём
		for (int j = 0; j < m; j++)
			table[i].a.push_back(i == j);
	}
}

// вывод таблицы
void SimplexMethod::Print() {
	int vars = variables.size();

	cout << endl;
	cout << "+-----+";

	for (int i = 0; i < vars; i++)
		cout << "-----------+";

	cout << endl;

	cout << "|  C  |";

	for (int i = 0; i < vars; i++)
		cout << " " << setw(9) << c[i] << " |";

	cout << endl;

	cout << "+-----+";

	for (int i = 0; i <= vars; i++)
		cout << "-----------+";

	cout << endl;

	cout << "|basis|";
	for (int i = 0; i < vars; i++)
		cout << "    x" << setw(2) << left << (i + 1) << "    |";

	cout << "     b     |" << endl;
	cout << "+-----+";

	for (int i = 0; i <= vars; i++)
		cout << "-----------+";

	cout << endl;

	for (int i = 0; i < m; i++) {
		cout << "| x" << setw(2) << left;

		if ((size_t)i < basis.size())
			cout << (basis[i] + 1);
		else
			cout << "?";

		cout << " |";

		for (size_t j = 0; j < table[i].a.size(); j++)
			cout << " " << setw(9) << (abs(table[i].a[j]) <= 1e-15 ? 0 : table[i].a[j]) << " |";

		cout << " " << setw(9) << (abs(table[i].b) <= 1e-15 ? 0 : table[i].b) << " |" << endl;
	}

	cout << "+-----+";

	for (int i = 0; i <= vars; i++)
		cout << "-----------+";

	cout << endl;

	if (!deltas.size())
		return;

	cout << "|  D  |";

	for (size_t i = 0; i < deltas.size(); i++)
		cout << " " << setw(9) << (abs(deltas[i]) <= 1e-15 ? 0 : deltas[i]) << " |";

	cout << endl;

	cout << "+-----+";

	for (int i = 0; i <= vars; i++)
		cout << "-----------+";

	cout << endl;


	cout << endl << endl << "X = (";

	for (int i = 0; i < n; i++) {
		bool inBasis = false;

		for (int j = 0; j < m; j++) {
			if (basis[j] == i) {
				cout << (abs(table[j].b) <= 1e-15 ? 0 : table[j].b) << " ";
				inBasis = true;
				solvation[i] = (abs(table[j].b) <= 1e-15 ? 0 : table[j].b);
				break;
			}
		}

		if (!inBasis) {
			solvation[i] = 0;
			cout << "0 ";
		}
	}

	cout << ")";

	cout << "\tF = " << deltas[n + m] << endl << endl;
}

// вычисление дельт
void SimplexMethod::CalculateDeltas() {
	deltas.clear(); // очищаем массив дельт

					// проходимся по всем переменным
	for (size_t i = 0; i <= variables.size(); i++) {
		double delta = 0;

		// вычилсяем дельту
		for (size_t j = 0; j < basis.size(); j++)
			delta += c[basis[j]] * (i < variables.size() ? table[j].a[i] : table[j].b);

		// вычитаем коэффициент функции
		if (i < variables.size())
			delta -= c[i];

		deltas.push_back(delta); // добавляем дельту в массив
	}
}

// вычисление номера минимальной дельты
int SimplexMethod::GetArgMaxDelta() {
	int imax = 0; // считаем, что первая дельта максимальна

				  // проходимся по всем дельтам
	for (size_t i = 1; i < deltas.size() - 1; i++)
		if (deltas[i] > deltas[imax]) // если дельта стала больше максимальной
			imax = i; // обновляем индекс максимума

	return imax; // возвращаем индекс максимума
}

// вычисление номера минимальной дельты
int SimplexMethod::GetArgMinDelta() {
	int imin = 0; // считаем, что первая дельта минимальная

				  // проходимся по всем дельтам
	for (size_t i = 1; i < deltas.size() - 1; i++)
		if (deltas[i] < deltas[imin]) // если дельта стала меньше минимальной
			imin = i; // обновляем индекс минимума

	return imin; // возвращаем индекс минимума
}

// решение ОЗЛП  max = 1 при минимизации max = -1 при максимизации
void SimplexMethod::Solve(int max) {
	while (true) {
		CalculateDeltas(); // рассчитываем дельты
		int jmax;

		// если минимизируем
		if (max == 1)
			jmax = GetArgMaxDelta(); // ищем индекс максимальной
									 // если максимизация
		else
			jmax = GetArgMinDelta(); // ищем индекс минимальной

		double maxDelta = deltas[jmax]; // получаем максимальную дельту

										// если она не положительна(или неотрицатльная для максимизации)
		if (maxDelta * max <= 0) {
			Print(); // выводим таблицу 
			isSolve = true;
			break; // и выходим
		}

		vector<double> Q(m); // создаём симплекс отношения
		int imin = -1;

		// идём по ограничениям
		for (int i = 0; i < m; i++) {
			if (table[i].a[jmax] == 0) { // если коэффициент равен 0
				Q[i] = 0; // то отношение равно нулю
			}
			else {
				Q[i] = table[i].b / table[i].a[jmax]; // вычисляем результат отношения

													  // если оно отрицательно, то идём дальше
				if (Q[i] < 0)
					continue;

				// иначе обновляем минимальное симплекс отношение
				if (imin == -1 || Q[i] < Q[imin])
					imin = i;
			}
		}

		if (imin == -1) {
			cout << "Solution does not exist" << endl;
			isSolve = false;
			break;
		}

		basis[imin] = jmax; // делаем переменную базисноц
		double pivot = table[imin].a[jmax]; // получаем опорный элемент

		MakeNewBasis(pivot, imin, jmax);
	}
}

// получение значения целевой функции по вектору
int SimplexMethod::GetF(vector<int> solve) {
	int res = 0;

	for (int i = 0; i < n; i++)
		res += c[i] * solve[i];

	return res;
}

// получение значения целевой функции по вектору
double SimplexMethod::GetF(vector<double> solve) {
	double res = 0;

	for (int i = 0; i < n; i++)
		res += c[i] * solve[i];

	return res;
}

// вывод вектора
void SimplexMethod::PrintVector(vector<int> vec) {
	cout << "(";

	for (size_t i = 0; i < vec.size() - 1; i++)
		cout << vec[i] << ",";

	cout << vec[vec.size() - 1] << ")";
}

// вывод вектора
void SimplexMethod::PrintVector(vector<double> vec) {
	cout << "(";

	for (size_t i = 0; i < vec.size() - 1; i++)
		cout << (abs(vec[i]) <= 1e-15 ? 0 : vec[i]) << ",";

	cout << vec[vec.size() - 1] << ")";
}

// проверка на то, что решение целое
int SimplexMethod::CheckForIntegerSolvation() {
	// если нет решения, то -1
	if (!isSolve)
		return -1;

	// проверка каждого элемента на то, что он целый
	for (int i = 0; i < n; i++)
		if (abs(solvation[i] - floor(solvation[i])) >= 1e-15)
			return 0; // если не целый, то возвращаем 0

	return 1; // если все целые, то возвращаем 1
}

// получение индекса максимального не целого числа 
int SimplexMethod::GetMaxFloatIndex() {
	double imax = 0; // индекс

	// находим индекс
	for (int i = 0; i < solvation.size(); i++) {
		if (solvation[i] - floor(solvation[i]) != 0 && solvation[i] > solvation[imax])
			imax = i;
	}

	return imax; // возвращаем
}

// получение условия
vector<double> SimplexMethod::GetRestriction(int index, double &res) {
	int j; // индекс строки, из которой будем брать условие

	// находим базисныю переменную по индексу
	for (int i = 0; i < basis.size(); i++) {
		if (i == index) {
			j = i;
			break;
		}
	}

	vector<double> restriction; // условие

	// записываем коэффициенты по правилу
	for (int i = 0; i < variables.size(); i++) {
		if (i < n)
			restriction.push_back(0);
		else
			restriction.push_back(-(table[j].a[i] - int(table[j].a[i])));
	}

	res = -(table[j].b - int(table[j].b)); 

	return restriction; // возвращаем условие
}

// поиск целочисленного решения иетодом Гомори
void SimplexMethod::IntegerSolve() {
	SimplexMethod restrictions; // таблица только для созранения условий задачи
	SimplexMethod methodGomori;

	// поиск решения(не целочисленного)
	cout << "Find solvation by simplex method" << endl;
	cout << "Initial table:" << endl;
	methodGomori.Print(); // вывод исходной таблицы

	cout << "Result table:" << endl;
	methodGomori.Solve(-1); // решение и вывод результата

	// пока решение есть 
	while (methodGomori.isSolve) {
		// если решение стало целочисленным, то выходим
		if (methodGomori.CheckForIntegerSolvation() == 1)
			break;

		int indexMaxFloat = methodGomori.GetMaxFloatIndex(); // находим максимальное не целое число среди решений

		cout << "Max real variable x" << (indexMaxFloat + 1) << " = " << methodGomori.solvation[indexMaxFloat] << endl;
		cout << "Build a cutting plane: ";

		double res;

		vector<double> restriction = methodGomori.GetRestriction(indexMaxFloat, res); // получаем ограничение

		cout << (restriction[0] >= 0 ? " " : "- ") << (abs(restriction[0]) <= 1e-15 ? 0 : abs(restriction[0])) << " ";

		// выводим ограничение
		for (int i = 1; i < restriction.size(); i++) {
			cout << (restriction[i] >= 0 ? "+ " : "- ") << (abs(restriction[i]) <= 1e-15 ? 0 : abs(restriction[i])) << " ";
		}

		cout << "= " << res << endl;

		restrictions.SetRestriction(restriction, res); // добавляем ограничение
		methodGomori = restrictions; // копируем данные

		cout << "Initial table:" << endl;
		methodGomori.Print(); // вывод исходной таблицы
		
		cout << "Result table:" << endl;
		methodGomori.Solve(-1); // решение и вывод результата	
	}

	// если нет решения, то выводим это
	if (!methodGomori.isSolve) {
		cout << "No integer solvation" << endl;
		return;
	}

	cout << "*********************************RESULT********************************************" << endl;

	cout << endl << "Solvation: ";
	PrintVector(methodGomori.solvation); // вывод решения
	cout << "\tF = " << methodGomori.GetF(methodGomori.solvation) << endl << endl; // и его функции цели
}

int main() {
	SimplexMethod method; // прямой метод
	method.IntegerSolve(); // запуск поиска
}
