#include <iostream> 
#include <cassert> 
#include <cmath> 
#include <iomanip> 

using namespace std;


// Класс для представления вектора в пространстве 
class vect {
public:
    int dim; // Размерность вектора 
    double* v; // Массив для хранения компонентов вектора 

    // Конструктор по умолчанию создает нулевой вектор в пространстве размерности dim 
    vect(int d) : dim(d), v(new double[d]) {}

    // Конструктор копирования создает новый вектор как копию другого вектора 
    vect(const vect& other) : dim(other.dim), v(new double[other.dim]) {
        // Копирование компонентов вектора 
        for (int i = 0; i < dim; ++i) {
            v[i] = other.v[i];
        }
    }

    // Деструктор освобождает память, выделенную для вектора 
    ~vect() {
        delete[] v;
    }

    // Оператор присваивания 
    vect& operator=(const vect& other) {
        // Проверка, присваиваем ли мы вектор самому себе 
        if (this != &other) {
            delete[] v; // Освобождение памяти, выделенной для текущего вектора 
            // Копирование размерности и выделение памяти для нового вектора 
            dim = other.dim;
            v = new double[dim];
            // Копирование компонентов вектора 
            for (int i = 0; i < dim; ++i) {
                v[i] = other.v[i];
            }
        }

        return *this; // Возвращаем ссылку на текущий вектор 
    }

    // Оператор сложения векторов 
    vect operator+(const vect& other) const {
        assert(dim == other.dim);  // Проверка, что размерности векторов совпадают 
        vect result(dim); // Создание результирующего вектора 
        // Покомпонентное сложение векторов 
        for (int i = 0; i < dim; ++i) {
            result.v[i] = v[i] + other.v[i];
        }
        return result;
    }

    // Оператор вычитания векторов 
    vect operator-(const vect& other) const {
        assert(dim == other.dim); // Проверка, что размерности векторов совпадают 
        vect result(dim);   // Создание результирующего вектора 
        // Покомпонентное вычитание векторов 
        for (int i = 0; i < dim; ++i) {
            result.v[i] = v[i] - other.v[i];
        }
        return result;
    }

    // Оператор унарного минуса 
    vect operator-() const {
        vect result(dim); // Создание результирующего вектора 

        // Покомпонентное умножение на -1 
        for (int i = 0; i < dim; ++i) {
            result.v[i] = -v[i];
        }

        return result;
    }

    // Оператор скалярного произведения векторов 
    double operator*(const vect& other) const {
        assert(dim == other.dim);  // Проверка, что размерности векторов совпадают 
        double result = 0;
        // Вычисление скалярного произведения 
        for (int i = 0; i < dim; ++i) {
            result += v[i] * other.v[i];
        }
        return result;
    }

    // Оператор умножения вектора на число 
    vect operator*(double k) const {
        vect result(dim);
        // Покомпонентное умножение на число 
        for (int i = 0; i < dim; ++i) {
            result.v[i] = v[i] * k;
        }
        return result;
    }
};

// Класс для представления матрицы 
class matr {
public:
    int dim; // Размерность матрицы 
    double** a;  // Массив для хранения элементов матрицы 

    // Конструктор по умолчанию создает нулевую матрицу размерности dim 
    matr(int d) : dim(d) {
        a = new double* [dim];  // Выделение памяти для массива указателей на строки матрицы 

        // Выделение памяти для каждой строки матрицы 
        for (int i = 0; i < dim; ++i) {
            a[i] = new double[dim];
        }
    }

    // Конструктор копирования создает новую матрицу как копию другой матрицы 
    matr(const matr& other) : dim(other.dim) {
        // Выделение памяти для массива указателей на строки матрицы 
        a = new double* [dim];
        // Выделение памяти для каждой строки матрицы и копирование элементов 
        for (int i = 0; i < dim; ++i) {
            a[i] = new double[dim];

       
        for (int j = 0; j < dim; ++j) {
            a[i][j] = other.a[i][j];
        }
        }
    }

    // Деструктор освобождает память, выделенную для матрицы 
    ~matr() {
        for (int i = 0; i < dim; ++i) {
            delete[] a[i];
        }
        delete[] a;
    }

    matr& operator=(const matr& other) {
        if (this != &other) {
            for (int i = 0; i < dim; ++i) {
                delete[] a[i];
            }
            delete[] a;
            dim = other.dim;
            a = new double* [dim];
            for (int i = 0; i < dim; ++i) {
                a[i] = new double[dim];
                for (int j = 0; j < dim; ++j) {
                    a[i][j] = other.a[i][j];
                }
            }
        }
        return *this;
    }

    matr operator+(const matr& other) const {
        assert(dim == other.dim);
        matr result(dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                result.a[i][j] = a[i][j] + other.a[i][j];
            }
        }
        return result;
    }

    matr operator-(const matr& other) const {
        assert(dim == other.dim);
        matr result(dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                result.a[i][j] = a[i][j] - other.a[i][j];
            }
        }
        return result;
    }

    matr operator-() const {
        matr result(dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                result.a[i][j] = -a[i][j];
            }
        }
        return result;
    }

    matr operator*(const matr& other) const {
        assert(dim == other.dim);
        matr result(dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                result.a[i][j] = 0;
                for (int k = 0; k < dim; ++k) {
                    result.a[i][j] += a[i][k] * other.a[k][j];
                }
            }
        }
        return result;
    }

    matr operator*(double k) const {
        matr result(dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                result.a[i][j] = a[i][j] * k;
            }
        }
        return result;
    }

    vect operator*(const vect& v) const {
        assert(dim == v.dim);
        vect result(dim);
        for (int i = 0; i < dim; ++i) {
            result.v[i] = 0;
            for (int j = 0; j < dim; ++j) {
                result.v[i] += a[i][j] * v.v[j];
            }
        }
        return result;
    }
};

//Оператор вывода вектора 
std::ostream& operator<<(std::ostream& os, const vect& v) {
    os << "vect(";
    for (int i = 0; i < v.dim; ++i) {
        os << v.v[i];
        if (i < v.dim - 1) {
            os << ", ";
        }
    }
    os << ")";
    return os;
}

//Оператор вывода матрицы 
std::ostream& operator<<(std::ostream& os, const matr& m) {
    os << "matr(";
    for (int i = 0; i < m.dim; ++i) {
        for (int j = 0; j < m.dim; ++j) {
            os << m.a[i][j];
            if (j < m.dim - 1) {
                os << ", ";
            }
        }
        if (i < m.dim - 1) {
            os << "; ";
        }
    }
    os << ")";
    return os;
}
// Функция для вычисления определителя матрицы A размера n x n
double determinant(const matr& A, int n) {

    double det = 0; // Инициализация переменной det - определитель
    int sign = 1; // Переменная sign для изменения знака при вычислении определителя

    for (int i = 0; i < n; i++) { // Цикл по столбцам матрицы
        matr submatrix(n); // Создание подматрицы размера (n-1) x (n-1)

        for (int j = 1; j < n; j++) { // Цикл по строкам подматрицы
            for (int k = 0; k < i; k++) {
                submatrix.a[j - 1][k] = A.a[j][k]; // Заполнение элементов подматрицы из матрицы A
            }
            for (int k = i + 1; k < n; k++) {
                submatrix.a[j - 1][k - 1] = A.a[j][k]; // Пропуск элемента, находящегося в выбранном столбце
            }
        }

        det += sign * A.a[0][i] * determinant(submatrix, n - 1); // Рекурсивный вызов функции для подматрицы
        sign = -sign; // Смена знака перед следующим слагаемым
    }

    return det; // Возвращаем вычисленное значение определителя
}

// Функция для вычисления обратной матрицы
matr inverse(matr A, int n) {
    double det = determinant(A, n); // Вычисление определителя матрицы A

    if (fabs(det) < 1e-10)
    {
        return false; // Возвращение false, если определитель близок к нулю, то есть матрица не имеет обратной
    }

    matr temp(n); // Создание временной матрицы temp для хранения результатов вычислений

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int sign = ((i + j) % 2 == 0) ? 1 : -1; // Определение знака для вычисления минора

            matr submatrix(n); // Создание подматрицы

            for (int k = 0; k < n; k++) {
                if (k != i) {
                    for (int l = 0; l < n; l++) {
                        if (l != j) {
                            submatrix.a[k - (k > i)][l - (l > j)] = A.a[k][l]; // Заполнение подматрицы элементами матрицы A за исключением i-ой строки и j-го столбца
                        }
                    }
                }
            }

            temp.a[j][i] = sign * determinant(submatrix, n - 1) / det; // Вычисление элемента обратной матрицы
        }
    }

    matr result(n); // Создание матрицы result для хранения результата обратной матрицы

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result.a[i][j] = temp.a[i][j]; // Копирование элементов из временной матрицы temp в результат обратной матрицы
        }
    }

    return result; // Возвращение обратной матрицы
}

//Функция exponet создает матрицу H, содержащую элементы, обратные элементам матрицы A
// Если элемент матрицы A равен нулю, соответствующий элемент матрицы H заполняется нулем. Возвращает матрицу H

matr exponet(matr A, int n) {
    matr H(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (A.a[i][j] == 0) H.a[i][j] = 0;
            H.a[i][j] = 1 / A.a[i][j];
        }
    }
    return H;
}

// Функция для решения системы линейных уравнений методом Якоби 
vect solveJacobi(const matr& A, const vect& B, const vect& X_0, int n) {
    matr H = exponet(A, n); //Сначала вычисляется матрица H, обратная матрице A, с использованием функции exponet
    matr I(n); //Создается единичная матрица I
    vect X(n);
    //Единичная матрица 
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) I.a[i][j] = 1;
            else I.a[i][j] = 0;
        }
    }
    //Затем в цикле производятся вычисления по методу Якоби: вычисляются две временные матрицы,
    // умножается матрица H на матрицу A, вычитается полученная матрица из единичной, 
    // вычисляются два вектора, их сумма присваивается X

    for (int i = 0; i < n; ++i) {
        matr m1 = H * A;
        matr m2 = I - m1;
        vect v1 = m2 * X_0;
        vect v2 = H * B;
        X = v1 + v2;
    }
    //Возвращается вектор X, содержащий решение системы уравнений
    return X;
}

int main() {

    setlocale(LC_ALL, "Russian"); //Устанавливает русскоязычную локаль для корректного отображения текста на кириллице
    int n; // Размерность системы 
    cout << "Введите размерность матрицы n = "; 
    cin >> n;

    matr A(n); //Объявление объекта матрицы A заданной размерности
    vect X_0(n); //Объявление объекта вектора X_0 заданной размерности

    srand(time(0)); //Генерация случайной матрицы A и вектора X_0 с помощью функции rand()

    // Генерация матрицы A 
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                A.a[i][j] = 100 * (rand()%50-20); // Увеличение диагональных элементов в 100 раз 
                if (A.a[i][j] == 0) A.a[i][j] = 100 * (rand() % 50-10);
            }
            else {
                A.a[i][j] = rand();
                if (A.a[i][j] == 0) A.a[i][j] = (rand() % 50-25);
            }

        }
    }
    //Вывод сгенерированной матрицы A
    cout << "Матрица A: ";
    cout << A << endl;
    //Заполнение вектора X_0 пользовательскими значениями
    for (int i = 0; i < n; ++i) {
        cout << "Введите X" << i << ": ";
        cin >> X_0.v[i];
    }
    //Вывод введенного пользователем вектора X_0
    cout << "Вектор X_0: ";
    cout << X_0 << endl;
    //Вычисление вектора B путем умножения матрицы A на вектор X_0
    vect B = A * X_0;
    cout << "Вектор B: ";
    cout << B << endl;

    // Вызов функции solveJacobi() для решения системы уравнений методом Якоби
    vect X = solveJacobi(A, B, X_0, n);

    // Вывод результатов 
    cout << "Решение системы X: ";
    cout << X << endl;

    return 0;
}