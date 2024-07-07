#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include <list>
#include "glut.h"

const double PI = acos(-1);

struct Point {
	//Fields
	double x;
	double y;

	//Constructors
	Point() = default;
	Point(double x, double y) {
		this->x = x;
		this->y = y;
	}
	Point(std::pair<double, double> other) {
		x = other.first;
		y = other.second;
	}

	//Methods
	double Dist_Between_Points(const Point& p) {
		return sqrt((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y));
	}

	//Operators
	void operator =(const Point& other) {
		this->x = other.x;
		this->y = other.y;
	}
	void operator =(const std::pair<double, double> other) {
		x = other.first;
		y = other.second;
	}

	//Bool operators
	bool operator ==(const Point& other) const {
		return (this->x == other.x) && (this->y == other.y);
	}
	bool operator !=(const Point& other) const {
		return (this->x != other.x) || (this->y != other.y);
	}
};
std::ostream& operator <<(std::ostream& cout, const Point& other) {
	std::cout << "X:" << other.x << " Y:" << other.y;
	return cout;
}
std::istream& operator >>(std::istream& cin, Point& other) {
	std::cout << "Введите координаты точки.\nX:";
	std::cin >> other.x;
	std::cout << "Y:";
	std::cin >> other.y;
	return cin;
}

class Line {
private:
	double _A;
	double _B;
	double _C;
public:
	//Constructors
	Line() = default;
	Line(const double& A, const double& B, const double& C) {
		this->_A = A;
		this->_B = B;
		this->_C = C;
	}
	//По двум точкам
	Line(const Point& p1, const Point& p2) {
		if (p1 == p2) throw std::overflow_error("Points must be different");
		this->_A = p2.y - p1.y;
		this->_B = p1.x - p2.x;
		this->_C = p2.x * p1.y - p1.x * p2.y;
	}
	//По угловому коэффициенту и сдвигу (y=kx+b)
	Line(const double& k, const double& b) {
		this->_A = -k;
		this->_B = 1;
		this->_C = -b;
	}
	//По точке и угловому коэффициенту
	Line(const Point& p, const double& k) {
		this->_A = -k;
		this->_B = 1;
		this->_C = -(_A * p.x + p.y);
	}

	//Getters
	double GetA() { return _A; }
	double GetB() { return _B; }
	double GetC() { return _C; }

	//Setters
	void Set_A(double new_A) { _A = new_A; }
	void Set_B(double new_B) { _B = new_B; }
	void Set_C(double new_C) { _C = new_C; }

	//Methods
	Point  point_of_intersection(Line& other) {
		Point ans;
		ans.x = (other.GetB() - _B + other.GetC() - _C) / (_A - other.GetA());
		ans.y = (other.GetA() - _A + other.GetC() - _C) / (_B - other.GetB());
		return ans;
	}

	//Operators
	void operator =(const Line& other) {
		this->_A = other._A;
		this->_B = other._B;
		this->_C = other._C;
	}

	//Bool operators
	bool operator ==(const Line& other) const {
		return (this->_A == other._A) && (this->_B == other._B) && (this->_C == other._C);
	}
	bool operator !=(const Line& other) const {
		return !(*this == other);
	}
};
std::istream& operator >>(std::istream& cin, Line& other) {
	double A, B, C;
	std::cout << "Введите коэффициенты из общего уравнения прямой (Ax+By+C=0).\nA:";
	std::cin >> A;
	std::cout << "B:";
	std::cin >> B;
	std::cout << "C:";
	std::cin >> C;
	other.Set_A(A); other.Set_B(B); other.Set_C(C);
	return cin;
}
std::ostream& operator <<(std::ostream& cout, Line& other) {
	std::cout << "A:" << other.GetA() << " B:" << other.GetB() << " C:" << other.GetC();
	return cout;
}

//-----------------------------------------------------------------------------------------------------------------------//

class Shape {
public:
	//Fields
	std::string Type = "Shape";

	//Getters
	std::string Get_Type() {
		return Type;
	}

	//glut
	virtual void Draw() = 0;
	virtual void rotate(Point center, double angle) = 0;
	virtual void reflex(Line axis) = 0;
	virtual void reflex(Point center) = 0;
	virtual void scale(Point center, double coef) = 0;

	//Math-functions
	virtual double  perimeter() = 0;
	virtual bool operator ==(Shape& other) = 0;
	virtual bool isCongruentTo(Shape& other) = 0;
	virtual double area() = 0;
	virtual bool containsPoint(Point& other) = 0;
	virtual bool isSimilarTo(Shape& other) = 0;
};

class Polygon : public Shape {
protected:
	std::vector<Point> _vertices;
	int _verticesCount;
	//std::string Type = "Polygon";
public:
	//Constructors
	Polygon() = default;
	Polygon(std::vector<Point> vertices) {
		_verticesCount = vertices.size();
		_vertices.resize(_verticesCount);
		for (int i = 0; i < _verticesCount; ++i) {
			_vertices[i] = vertices[i];
		}
		Type = "Polygon";
	}
	Polygon(std::initializer_list<Point> list) {
		for (auto& x : list) {
			_vertices.push_back(x);
		}
		_verticesCount = _vertices.size();
		Type = "Polygon";
	}

	//Getters
	int Get_verticesCount() {
		return _verticesCount;
	}
	std::vector<Point>& Get_vertices() {
		return _vertices;
	}

	//Methods
	bool isConvex() {
		if (_verticesCount == 3) return true;
		for (int i = 0; i < _verticesCount - 1; ++i) {
			Line wr_line(_vertices[i], _vertices[i + 1]);
			for (int j = i + 2; j < _verticesCount - 1; ++j) {
				double fx1; // Variable to store a * x1 + b * y1 - c
				double fx2; // Variable to store a * x2 + b * y2 - c

				fx1 = wr_line.GetA() * _vertices[j].x + wr_line.GetB() * _vertices[j].y - wr_line.GetC();
				fx2 = wr_line.GetA() * _vertices[j + 1].x + wr_line.GetB() * _vertices[j + 1].y - wr_line.GetC();

				if (fx1 * fx2 < 0) return false;
			}
		}
		return true;
	}

	//Override glut-methods
	void Draw() override {
		glLineWidth(2);
		glColor3ub(0, 0, 0);
		for (int i = 0; i < _verticesCount - 1; ++i) {
			glBegin(GL_LINES);
			glVertex2d(_vertices[i].x, _vertices[i].y);
			glVertex2d(_vertices[i + 1].x, _vertices[i + 1].y);
			glEnd();
		}
		glBegin(GL_LINES);
		glVertex2d(_vertices[0].x, _vertices[0].y);
		glVertex2d(_vertices[_verticesCount - 1].x, _vertices[_verticesCount - 1].y);
		glEnd();
	}
	void rotate(Point center, double angle) override {
		std::vector<Point> v(_vertices.size());
		for (int i = 0; i < _vertices.size(); ++i) {
			v[i].x = center.x + (_vertices[i].x - center.x) * cos(angle / 180 * PI) - (_vertices[i].y - center.y) * sin(angle / 180 * PI);
			v[i].y = center.y + (_vertices[i].x - center.x) * sin(angle / 180 * PI) + (_vertices[i].y - center.y) * cos(angle / 180 * PI);
		}
		_vertices = v;
	}
	void reflex(Line axis) override {
		int sz = _vertices.size();
		std::vector<Point> v = _vertices;
		if (axis.GetA() == 0) {
			for (int i = 0; i < sz; ++i) {
				_vertices[i].y = 2 * (-axis.GetC() / axis.GetB()) - _vertices[i].y;
			}
		}
		else if (axis.GetB() == 0) {
			for (int i = 0; i < sz; ++i) {
				_vertices[i].x = 2 * (-axis.GetC() / axis.GetB()) - _vertices[i].x;
			}
		}
		else {
			double k = -axis.GetA() / axis.GetB(), c = -axis.GetC() / axis.GetB();
			for (int i = 0; i < sz; ++i) {
				v[i].x = 2 * (_vertices[i].x - k * (c - _vertices[i].y)) / (1 + k * k) - _vertices[i].x;
				v[i].y = (k * (2 * _vertices[i].x + k * _vertices[i].y) + 2 * c - _vertices[i].y) / (1 + k * k);
			}
			_vertices = v;
		}
	}
	void reflex(Point center) override {
		for (int i = 0; i < _vertices.size(); ++i) {
			_vertices[i].x += 2 * (center.x - _vertices[i].x);
			_vertices[i].y += 2 * (center.y - _vertices[i].y);
		}
	}
	void scale(Point center, double coef) override {
		for (int i = 0; i < _vertices.size(); ++i) {
			_vertices[i].x = coef * (_vertices[i].x - center.x) + center.x;
			_vertices[i].y = coef * (_vertices[i].y - center.y) + center.y;
		}
	}

	//Override math-methods
	double perimeter() override {
		double ans = 0;
		for (int i = 0; i < _verticesCount - 1; ++i) {
			ans += _vertices[i].Dist_Between_Points(_vertices[i + 1]);
		}
		ans += _vertices[0].Dist_Between_Points(_vertices[_verticesCount - 1]);
		return ans;
	}
	bool operator ==(Shape& other) override {
		if (other.Type == "Ellipse" || other.Type == "Circle") {
			return false;
		}
		Polygon& wr = dynamic_cast<Polygon&>(other);
		if (_verticesCount != wr._verticesCount) return false;
		for (int i = 0; i < _verticesCount; ++i) {
			if (_vertices[i] != wr._vertices[i]) return false;
		}
		return true;
	}
	bool isCongruentTo(Shape& other) override {
		if (other.Type == "Ellipse" || other.Type == "Circle") {
			return false;
		}
		Polygon& wr = dynamic_cast<Polygon&>(other);
		if (_verticesCount != wr._verticesCount) return false;

		std::vector<double> arr1, arr2;
		for (int i = 0; i < _verticesCount - 1; ++i) {
			arr1.push_back(_vertices[i].Dist_Between_Points(_vertices[i + 1]));
			arr2.push_back(wr._vertices[i].Dist_Between_Points(wr._vertices[i + 1]));
		}
		arr1.push_back(_vertices[0].Dist_Between_Points(_vertices[_verticesCount - 1]));
		arr2.push_back(wr._vertices[0].Dist_Between_Points(wr._vertices[wr._verticesCount - 1]));
		std::sort(arr1.begin(), arr1.end());
		std::sort(arr2.begin(), arr2.end());
		for (int i = 0; i < _verticesCount; ++i) {
			if (arr1[i] != arr2[i]) return false;
		}
		return true;
	}
	double area() override {
		int sz = _verticesCount;
		double s = 0;
		for (int i = 0; i < sz - 1; ++i) {
			s += _vertices[i].x * _vertices[i + 1].y;
		}
		s += _vertices[sz - 1].x * _vertices[0].y;
		for (int i = 0; i < sz - 1; ++i) {
			s -= _vertices[i + 1].x * _vertices[i].y;
		}
		s -= _vertices[0].x * _vertices[sz - 1].y;
		s = abs(s);
		s /= 2;
		return s;
	}
	bool containsPoint(Point& other) override {
		//Если луч, выпущенный из точки пересекает четное количество сторон, то
		//точка не лежит внутри многоугольника. Иначе точка лежит в многоугольнике
		bool res = false;
		int j = _verticesCount - 1;
		for (int i = 0; i < _vertices.size(); i++) {
			if ((_vertices[i].y < other.y && _vertices[j].y >= other.y || _vertices[j].y < other.y && _vertices[i].y >= other.y) &&
				(_vertices[i].x + (other.y - _vertices[i].y) / (_vertices[j].y - _vertices[i].y) * (_vertices[j].x - _vertices[i].x) < other.x))
				res = !res;
			j = i;
		}
		return res;
	}
	bool isSimilarTo(Shape& other) override {
		if (other.Type == "Ellipse" || other.Type == "Circle") return false;

		Polygon& a = dynamic_cast<Polygon&>(other);

		if (_vertices.size() != a._vertices.size()) return false;
		if (_vertices.size() == 1) return true;
		int sz = _verticesCount;

		double x1 = _vertices[0].x - _vertices[sz - 1].x;
		double y1 = _vertices[0].y - _vertices[sz - 1].y;
		double x2 = _vertices[1].x - _vertices[0].x;
		double y2 = _vertices[1].y - _vertices[0].y;
		std::vector<double> vector_of_skalyars(sz);
		std::vector<double> vector_of_lengths(sz);
		vector_of_skalyars[0] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
		vector_of_lengths[0] = sqrt(x1 * x1 + y1 * y1);
		for (int i = 2; i < sz; i++) {
			x1 = x2;
			y1 = y2;
			x2 = _vertices[i].x - _vertices[i - 1].x;
			y2 = _vertices[i].y - _vertices[i - 1].y;
			vector_of_skalyars[i - 1] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
			vector_of_lengths[i - 1] = sqrt(x1 * x1 + y1 * y1);
		}

		x1 = x2;
		y1 = y2;
		x2 = _vertices[0].x - _vertices[sz - 1].x;
		y2 = _vertices[0].y - _vertices[sz - 1].y;
		vector_of_skalyars[sz - 1] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
		vector_of_lengths[sz - 1] = sqrt(x1 * x1 + y1 * y1);

		x1 = a._vertices[0].x - a._vertices[sz - 1].x;
		y1 = a._vertices[0].y - a._vertices[sz - 1].y;
		x2 = a._vertices[1].x - a._vertices[0].x;
		y2 = a._vertices[1].y - a._vertices[0].y;
		std::vector<double> vector_of_skalyars_sec(sz);
		std::vector<double> vector_of_lengths_sec(sz);
		vector_of_skalyars_sec[0] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
		vector_of_lengths_sec[0] = sqrt(x1 * x1 + y1 * y1);
		for (int i = 2; i < sz; i++) {
			x1 = x2;
			y1 = y2;
			x2 = a._vertices[i].x - a._vertices[i - 1].x;
			y2 = a._vertices[i].y - a._vertices[i - 1].y;
			vector_of_skalyars_sec[i - 1] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
			vector_of_lengths_sec[i - 1] = sqrt(x1 * x1 + y1 * y1);
		}

		x1 = x2;
		y1 = y2;
		x2 = a._vertices[0].x - a._vertices[sz - 1].x;
		y2 = a._vertices[0].y - a._vertices[sz - 1].y;
		vector_of_skalyars_sec[sz - 1] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
		vector_of_lengths_sec[sz - 1] = sqrt(x1 * x1 + y1 * y1);


		bool flag = true, id;
		double k;
		for (int i = 0; i < sz; i++) {
			flag = false;
			if (vector_of_skalyars_sec[0] == vector_of_skalyars[i]) {
				flag = true;
				id = i;
				k = vector_of_lengths_sec[0] / vector_of_lengths[i];
				for (int j = i + 1, cnt = 1; cnt < sz; ++cnt, ++j) {
					if (j == sz) j = 0;
					if (vector_of_skalyars_sec[cnt] != vector_of_skalyars[j] || k != (vector_of_lengths_sec[cnt] / vector_of_lengths[j])) {
						flag = false;
						break;
					}
				}
			}
			if (flag) break;
		}

		return flag;
	}
};
std::ostream& operator <<(std::ostream& cout, Polygon& other) {
	std::cout << "Координаты вершин многоугольника в порядке обхода:\n";
	for (int i = 0; i < other.Get_verticesCount(); ++i) {
		std::cout << i + 1 << ") " << other.Get_vertices()[i] << '\n';
	}
	return cout;
}

class Ellipse : public Shape {
protected:
	std::pair<Point, Point> _focuses;
	double _a, _b, _c;
public:
	//Constructors
	Ellipse() = default;
	Ellipse(const Point& f1, const Point& f2, const double& dist) {
		_focuses.first = f1;
		_focuses.second = f2;
		_a = dist / 2;
		_c = sqrt(pow(_focuses.first.x - _focuses.second.x, 2) + pow(_focuses.first.y - _focuses.second.y, 2)) / 2;
		_b = sqrt(_a * _a - _c * _c);
		Type = "Ellipse";
	}

	//Getters
	std::pair<Point, Point> Get_focuses() {
		return _focuses;
	}
	double Get_sumOfDistanse() {
		return _a * 2.;
	}
	std::string Get_Type() {
		return Type;
	}

	//Operators
	void operator =(const Ellipse& other) {
		_focuses = other._focuses;
		_a = other._a;
		_b = other._b;
		_c = other._c;
	}

	//Methods
	Point centre() {
		Point ans;
		ans.x = (_focuses.first.x + _focuses.second.x) / 2;
		ans.y = (_focuses.first.y + _focuses.second.y) / 2;
		return ans;
	}
	double eccentrisity() {
		return _c / _a;
	}
	std::pair<Line, Line> directrices() {
		return std::make_pair(Line(0, _b / (_c / _a)), Line(0, -_b / (_c / _a)));
	}

	//Override glut-methods
	void Draw() override {
		Point f1 = _focuses.first, f2 = _focuses.second;
		double Center_X = (f1.x + f2.x) / 2;
		double Center_Y = (f1.y + f2.y) / 2;
		const double step = float(2 * PI / 360);
		double r_x = _a, r_y = _b;
		Line l(f1, f2);
		double k = -l.GetA() / l.GetB();
		glLineWidth(2);
		glBegin(GL_LINE_STRIP);
		glColor3f(0, 0, 0);
		double x, y, x1, y1;
		for (double angle = 0; angle < double(2 * PI); angle += step) {
			x = r_x * cosf(angle) + Center_X;
			y = r_y * sinf(angle) + Center_Y;
			if (k != 0) {
				x1 = Center_X + (x - Center_X) * cos(atan(k)) - (y - Center_Y) * sin(atan(k));
				y1 = Center_Y + (x - Center_X) * sin(atan(k)) + (y - Center_Y) * cos(atan(k));
			}
			else {
				x1 = x;
				y1 = y;
			}
			glVertex2f(x1, y1);
		}
		glEnd();
		glFlush();
	}
	void rotate(Point center, double angle) override {
		Point g1, g2;
		Point f1 = _focuses.first, f2 = _focuses.second;
		g1.x = center.x + (f1.x - center.x) * cos(angle / 180 * PI) - (f1.y - center.y) * sin(angle / 180 * PI);
		g1.y = center.y + (f1.x - center.x) * sin(angle / 180 * PI) + (f1.y - center.y) * cos(angle / 180 * PI);

		g2.x = center.x + (f2.x - center.x) * cos(angle / 180 * PI) - (f2.y - center.y) * sin(angle / 180 * PI);
		g2.y = center.y + (f2.x - center.x) * sin(angle / 180 * PI) + (f2.y - center.y) * cos(angle / 180 * PI);

		_focuses.first = g1;
		_focuses.second = g2;
	}
	void reflex(Point center) override {
		Point f1 = _focuses.first, f2 = _focuses.second;

		_focuses.first.x += 2 * (center.x - f1.x);
		_focuses.first.y += 2 * (center.y - f1.y);

		_focuses.second.x += 2 * (center.x - f2.x);
		_focuses.second.y += 2 * (center.y - f2.y);
	}
	void reflex(Line axis) override {
		Point f1 = _focuses.first, f2 = _focuses.second;
		if (axis.GetA() == 0) {
			_focuses.first.y = 2 * (-axis.GetC() / axis.GetB()) - f1.y;
			_focuses.second.y = 2 * (-axis.GetC() / axis.GetB()) - f2.y;
		}
		else if (axis.GetB() == 0) {
			_focuses.first.x = 2 * (-axis.GetC() / axis.GetA()) - f1.x;
			_focuses.second.x = 2 * (-axis.GetC() / axis.GetA()) - f2.x;
		}
		else {
			double k = -axis.GetA() / axis.GetB(), c = -axis.GetC() / axis.GetB();
			Point g1, g2;
			g1.x = 2 * (f1.x - k * (c - f1.y)) / (1 + k * k) - f1.x;
			g1.y = (k * (2 * f1.x + k * f1.y) + 2 * c - f1.y) / (1 + k * k);
			g2.x = 2 * (f2.x - k * (c - f2.y)) / (1 + k * k) - f2.x;
			g2.y = (k * (2 * f2.x + k * f2.y) + 2 * c - f2.y) / (1 + k * k);
			_focuses.first = g1;
			_focuses.second = g2;
		}
	}
	void scale(Point center, double coef) override {
		Point f1 = _focuses.first, f2 = _focuses.second;
		f1.x = coef * (f1.x - center.x) + center.x;
		f1.y = coef * (f1.y - center.y) + center.y;
		f2.x = coef * (f2.x - center.x) + center.x;
		f2.y = coef * (f2.y - center.y) + center.y;
		_a *= coef;
		_b *= coef;
		_c *= coef;
	}

	//Override math-methods
	double perimeter() override {
		return 2 * PI * sqrt((_a * _a + _b * _b) / 2);
	}
	bool operator ==(Shape& other) override {
		if (other.Type != "Ellipse" && other.Type != "Circle") return false;
		Ellipse& wr = dynamic_cast<Ellipse&>(other);
		if (_focuses != wr._focuses || _a != wr._a || _b != wr._b || _c != wr._c) return false;
		return true;
	}
	bool isCongruentTo(Shape& other) override {
		if (other.Type != "Ellipse" && other.Type != "Circle") return false;
		Ellipse& wr = dynamic_cast<Ellipse&>(other);
		if (this->perimeter() != wr.perimeter()) return false;
		if (_a != wr._a || _b != wr._b || _c != wr._c) return false;
		return true;
	}
	double area() override {
		return PI * _a * _b;
	}
	bool containsPoint(Point& other) override {
		double r1, r2;
		r1 = sqrt(pow(_focuses.first.x - other.x, 2) + pow(_focuses.first.y - other.y, 2));
		r2 = sqrt(pow(_focuses.second.x - other.x, 2) + pow(_focuses.second.y - other.y, 2));
		return (r1 + r2 <= 2 * _a);
	}
	bool isSimilarTo(Shape& other) override {
		if (other.Type != "Ellipse" && other.Type != "Circle") return false;
		Ellipse& an = dynamic_cast<Ellipse&>(other);
		return (_a / an._a == _b / an._b && _a / an._a == _c / an._c);
	}
};
std::ostream& operator <<(std::ostream& cout, Ellipse& other) {
	std::cout << "Фокусы элипса и его сумма расстояний от фокусов до контура:\nF1) " << other.Get_focuses().first << "\nF2) " << other.Get_focuses().second << "\nDistanse) " << other.Get_sumOfDistanse();
	return cout;
}

class Circle : public Ellipse {
public:
	//Constructors
	Circle() : Ellipse() {
		Type = "Circle";
	}
	Circle(const Point& centre, double r) : Ellipse(centre, centre, r * 2.) {
		Type = "Circle";
	}

	//Operators
	void operator =(const Circle& other) {
		this->_a = other._a;
		this->_focuses = other._focuses;
		this->_b = other._b;
		this->_c = other._c;
	}

	//Getters
	Point centre() {
		return _focuses.first;
	}
	double radius() {
		return _a;
	}
	std::string Get_Type() {
		return Type;
	}

	//Glut-methods
	void Draw() {
		int pointCount = 360;
		const float step = float(2 * PI / pointCount);
		glLineWidth(2);
		glBegin(GL_LINE_LOOP);
		glColor3ub(0,0,0);
		for (double angle = 0; angle < 2. * PI; angle += step) {
			double dx = _a * cosf(angle);
			double dy = _a * sinf(angle);
			glVertex2f(_focuses.first.x + dx, _focuses.first.y + dy);
		}
		glEnd();
		glFlush();
	}
};
std::ostream& operator <<(std::ostream& cout, Circle& other) {
	std::cout << "Координаты центра окружности и её радиус:\nЦентр) " << other.Get_focuses().first << "\nРадиус) " << other.radius();
	return cout;
}

class Rectangle : public Polygon {
public:
	//Constructors
	Rectangle() = default;
	Rectangle(std::vector<Point> vertices) : Polygon(vertices) {
		Type = "Rectangle";
	}
	Rectangle(std::initializer_list<Point> lst) : Polygon(lst) {
		Type = "Rectangle";
	}
	Rectangle(Point p1, Point p2, double coef) {
		Type = "Rectangle";

		if (coef < 1) {
			coef = 1 / coef;
		}
		Point c((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
		Polygon pn({ p1, p2 });
		pn.rotate(c, 180 / (coef + 1) * coef);
		std::vector<Point> ver = pn.Get_vertices();
		_vertices.resize(4);
		_vertices[0] = p1;
		_vertices[1] = ver[0];
		_vertices[2] = p2;
		_vertices[3] = ver[1];
		if (_vertices[0].Dist_Between_Points(_vertices[3]) > _vertices[3].Dist_Between_Points(_vertices[2])) {
			this->reflex(Line(p1, p2));
		}
		_verticesCount = 4;
	}

	//Getters
	std::string Get_Type() {
		return Type;
	}

	//Methods
	std::pair<Line, Line> diagonals() {
		Line l1(_vertices[0], _vertices[2]);
		Line l2(_vertices[1], _vertices[3]);
		std::pair<Line, Line> ans = { l1, l2 };
		return ans;
	}
	Point centre() {
		Rectangle copy = *this;
		std::pair<Line, Line> diagonals = copy.diagonals();
		return diagonals.first.point_of_intersection(diagonals.second);
	}
	double area() {
		double ans = _vertices[0].Dist_Between_Points(_vertices[1]) * _vertices[1].Dist_Between_Points(_vertices[2]);
		return ans;
	}
};

class Square : public Rectangle {
public:
	//Constructors
	Square() : Rectangle() {
		Type = "Square";
	}
	Square(std::vector<Point> arr) : Rectangle(arr) {
		Type = "Square";
	}
	Square(const Point& p1, const Point& p2) : Rectangle(p1, p2, 1) {
		Type = "Square";
	}
	Square(std::initializer_list<Point> lst) : Rectangle(lst) {
		Type = "Square";
	}

	//Getters
	std::string Get_Type() {
		return Type;
	}

	//Methods
	Circle inscribedCircle() {
		double wr_r = _vertices[0].Dist_Between_Points(_vertices[1]) / 2.;
		Point wr_p((_vertices[0].x + _vertices[2].x) / 2, (_vertices[0].y + _vertices[2].y) / 2);
		Circle ans(wr_p, wr_r);
		return ans;
	}
	Circle circumscribedCircle() {
		Point wr_p((_vertices[0].x + _vertices[2].x) / 2, (_vertices[0].y + _vertices[2].y) / 2);
		double wr_r = _vertices[0].Dist_Between_Points(_vertices[2]) / 2;
		Circle ans(wr_p, wr_r);
		return ans;
	}
};

class Triangle : public Polygon {
public:
	//Constructors
	Triangle() : Polygon() {
		Type = "Triangle";
	}
	Triangle(std::vector<Point> arr) : Polygon(arr) {
		Type = "Triangle";
	}
	Triangle(std::initializer_list<Point> lst) : Polygon(lst) {
		Type = "Triangle";
	}

	//Getters
	std::string Get_Type() {
		return Type;
	}

	//Methods
	Circle circumscribedCircle() {
		Point A = _vertices[0], B = _vertices[1], C = _vertices[2];
		double D = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));

		double centerX_of_Circle = ((pow(A.x, 2) + pow(A.y, 2)) * (B.y - C.y) +
			(pow(B.x, 2) + pow(B.y, 2)) * (C.y - A.y) +
			(pow(C.x, 2) + pow(C.y, 2)) * (A.y - B.y)) / D;

		double centerY_of_Circle = ((pow(A.x, 2) + pow(A.y, 2)) * (C.x - B.x) +
			(pow(B.x, 2) + pow(B.y, 2)) * (A.x - C.x) +
			(pow(C.x, 2) + pow(C.y, 2)) * (B.x - A.x)) / D;

		Point p(centerX_of_Circle, centerY_of_Circle);
		double r = A.Dist_Between_Points(p);
		return Circle(p, r);
	}
	Circle inscribedCircle() {
		Point A = _vertices[0], B = _vertices[1], C = _vertices[2];
		double a = B.Dist_Between_Points(C);
		double b = A.Dist_Between_Points(C);
		double c = A.Dist_Between_Points(B);
		double wr = a + b + c;

		double centerX_of_Circle = (a * A.x + b * B.x + c * C.x) / wr;
		double centerY_of_Circle = (a * A.y + b * B.y + c * C.y) / wr;
		Point p(centerX_of_Circle, centerY_of_Circle);

		double s = this->area();
		double r = s / ((a + b + c) / 2);
		return Circle(p, r);
	}
	Point centroid() {//Точка пересечения медиан
		Point wr_point((_vertices[0].x + _vertices[1].x) / 2, (_vertices[0].y + _vertices[1].y) / 2);
		return Point((_vertices[2].x + 2 * wr_point.x) / 3, (_vertices[2].y + 2 * wr_point.y) / 3);
	}
	Point orthocenter() {//Точка пересечения высот
		Point ans;
		Point A = _vertices[0], B = _vertices[1], C = _vertices[2];
		double a = B.Dist_Between_Points(C);
		double b = A.Dist_Between_Points(C);
		double c = A.Dist_Between_Points(B);
		ans.x = (A.x * a * a + B.x * b * b + C.x * c * c) / (a * a + b * b + c * c);
		ans.y = (A.y * a * a + B.y * b * b + C.y * c * c) / (a * a + b * b + c * c);
		return ans;
	}
	Line EulerLine() {
		return Line(this->centroid(), this->orthocenter());
	}
	Circle ninePointsCircle() {
		Point wr1 = (this->inscribedCircle()).centre();
		Point wr2 = this->orthocenter();
		Point wr3((wr1.x + wr2.x) / 2, (wr1.y + wr2.y) / 2);
		Point wr4((_vertices[0].x + _vertices[1].x) / 2, (_vertices[0].y + _vertices[1].y) / 2);
		double r = wr3.Dist_Between_Points(wr4);
		return Circle(wr3, r);
	}
};

//-----------------------------------------------------------------------------------------------------------------------//

int Max_Coordinate = 10;
int Flag_Conditions = 1;
//1: NEW
//2: ADD
//3: CHA -> 4: MAT -> 5: ==


struct Buttom {
	std::string name;
	double x = -Max_Coordinate + .5, y = Max_Coordinate - .5;
	Polygon body;
	double width = 2, hight = 1;
	double R = 204, G = 204, B = 255;

	Buttom(std::string name, double bord) {
		this->name = name;
		y -= bord;
		Point p1(x,y);
		Point p2(x + width, y);
		Point p3(x + width, y - hight);
		Point p4(x, y - hight);
		body = { p1, p2, p3, p4 };
	}
	Buttom(std::string name, double bord, double R, double G, double B) {
		this->name = name;
		y -= bord;
		this->R = R;
		this->G = B;
		this->B = B;
		Point p1(x, y);
		Point p2(x + width, y);
		Point p3(x + width, y - hight);
		Point p4(x, y - hight);
		body = { p1, p2, p3, p4 };
	}

	void Draw() {
		glColor3ub(R, G, B);
		glBegin(GL_QUADS);
		glVertex2d(x, y);
		glVertex2d(x + width, y);
		glVertex2d(x + width, y - hight);
		glVertex2d(x, y - hight);
		glEnd();

		glLineWidth(1);
		glColor3ub(0, 0, 0);
		glBegin(GL_LINES);
		glVertex2d(x, y);
		glVertex2d(x + width, y);

		glVertex2d(x + width, y);
		glVertex2d(x + width, y - hight);

		glVertex2d(x + width, y - hight);
		glVertex2d(x, y - hight);

		glVertex2d(x, y - hight);
		glVertex2d(x, y);
		glEnd();

		glColor3ub(0, 0, 0);
		glLineWidth(3);
		if (name == "Esc") {
			glBegin(GL_LINES);
			glVertex2d(x + .5, y - .5);
			glVertex2d(x + 1.5, y - .5);

			glVertex2d(x + .5, y - .5);
			glVertex2d(x + .7, y - .3);

			glVertex2d(x + .5, y - .5);
			glVertex2d(x + .7, y - .7);
			glEnd();
			return;
		}

		if (name == "New") {
			glRasterPos2f(x + .5, y - .7);
			for (int i = 0; i < name.size(); ++i) {
				glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, name[i]);
			}
			glFlush();
			return;
		}

		glRasterPos2f(x + .01, y - .7);
		for (int i = 0; i < name.size(); ++i) {
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, name[i]);
		}
		glFlush();
		return;
	}
};

std::vector<Shape*> Shapes;
std::vector<Buttom> Buttoms;
std::vector<Buttom> Buttoms_Func;
std::vector<Buttom> Buttoms_Math_Func;

Shape* Selected_Shape;
Shape* Selected_Shape_Second;
std::string flag_name;

void Draw_Coordinate_Axes() {
	glBegin(GL_LINES);
	glVertex2d(-.1, 1.);
	glVertex2d(.1, 1);

	glVertex2d(1., -.1);
	glVertex2d(1, .1);

	glVertex2d(-1., -.1);
	glVertex2d(-1, .1);

	glVertex2d(-.1, -1.);
	glVertex2d(.1, -1.);
	glEnd();

	glColor3ub(100, 100, 100);
	glLineWidth(1);
	glBegin(GL_LINES);
	for (double i = -10.0; i <= 10.0; i += 1.0) {
		glVertex2d(i, -10.0); // Линии по оси X
		glVertex2d(i, 10.0);
		glVertex2d(-10.0, i); // Линии по оси Y
		glVertex2d(10.0, i);
	}
	glEnd();

	glLineWidth(1);
	glColor3ub(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2d(0., -1000000.);
	glVertex2d(0., 1000000);

	glVertex2d(-1000000., 0.);
	glVertex2d(1000000, 0.);
	glEnd();
}

void Create_Shape(std::string name) {
	if (name == "Polygon") {
		std::vector<Point> vec;
		std::cout << "Введите количество вершин: ";
		int N;
		std::cin >> N;
		while (N <= 2) {
			std::cout << "Некорректный ввод. Попробуйте снова: ";
			std::cin >> N;
		}
		std::cout << "Введите координаты вершин в порядке обхода: " << '\n';
		for (int i = 0; i < N; ++i) {
			Point x;
			std::cin >> x;
			vec.push_back(x);
		}
		std::cout << "Фигура отрисована" << std::endl << std::endl;
		Polygon* x = new Polygon(vec);
		Shapes.push_back(x);
		return;
	}
	if (name == "Ellipse") {
		Point f1, f2;
		std::cout << "Задайте фокусы эллипса: " << '\n';
		std::cin >> f1 >> f2;
		double t;
		std::cout << "Задайте сумму расстояний от фокусов до кривой: ";
		std::cin >> t;
		std::cout << "Фигура отрисована\nИнформация о фигуре:" << std::endl;
		Ellipse* x = new Ellipse(f1, f2, t);
		Shapes.push_back(x);
		std::cout << "Perimetr: " << x->perimeter() << '\n' << "Area: " << x->area() << std::endl << std::endl;
		return;
	}
	if (name == "Circle") {
		Point p;
		std::cout << "Задайте центр окружности: ";
		std::cin >> p;
		double r;
		std::cout << "Введите радиус окружности: ";
		std::cin >> r;
		std::cout << "Фигура отрисована" << std::endl << std::endl;
		Circle* x = new Circle(p, r);
		Shapes.push_back(x);
		return;
	}
	if (name == "Rectangle") {
		std::vector<Point> vec;
		std::cout << "Введите координаты вершин в порядке обхода: " << '\n';
		for (int i = 0; i < 4; ++i) {
			Point x;
			std::cin >> x;
			vec.push_back(x);
		}
		std::cout << "Фигура отрисована" << std::endl << std::endl;
		Rectangle* x = new Rectangle(vec);
		Shapes.push_back(x);
		return;
	}
	if (name == "Square") {
		std::vector<Point> vec;
		std::cout << "Введите координаты вершин в порядке обхода: " << '\n';
		for (int i = 0; i < 4; ++i) {
			Point x;
			std::cin >> x;
			vec.push_back(x);
		}
		std::cout << "Фигура отрисована" << std::endl << std::endl;
		Square* x = new Square(vec);
		Shapes.push_back(x);
		return;
	}
	if (name == "Triangle") {
		std::vector<Point> vec;
		std::cout << "Введите координаты вершин в порядке обхода: " << '\n';
		for (int i = 0; i < 3; ++i) {
			Point x;
			std::cin >> x;
			vec.push_back(x);
		}
		std::cout << "Фигура отрисована" << std::endl << std::endl;
		Triangle* x = new Triangle(vec);
		Shapes.push_back(x);
		return;
	}
}
void ReNew_Information(Buttom bt) {
	int indx = 0;
	for (int i = 0; i < Shapes.size(); ++i) {
		if (Selected_Shape == Shapes[i]) indx = i;
	}
	if (bt.name == "Del") {
		Shapes.erase(Shapes.begin() + indx);
		Flag_Conditions = 1;
		std::cout << "Удаление выполнено\n\n";
		return;
	}
	if (bt.name == "Rotate") {
		double angle;
		Point p;
		std::cout << "Введите данные для функции Rotate\nВведите угол поворота: ";
		std::cin >> angle;
		std::cin >> p;

		Selected_Shape->rotate(p, angle);
		Flag_Conditions = 1;
		std::cout << "Изменения внесены\n\n";
		return;
	}
	if (bt.name == "Scale") {
		Point p;
		double k;
		std::cout << "Введите данные для гомотетии\nВведите коэффициент для гомотетии: ";
		std::cin >> k;
		std::cin >> p;
		Selected_Shape->scale(p, k);
		Flag_Conditions = 1;
		std::cout << "Изменения внесены\n\n";
		return;
	}
	if (bt.name == "P_Reflex") {
		Point p;
		std::cout << "Введите данные для функции Rexlex относительно точки\n";
		std::cin >> p;
		Selected_Shape->reflex(p);
		Flag_Conditions = 1;
		std::cout << "Изменения внесены\n\n";
		return;
	}
	if (bt.name == "L_Reflex") {
		
		std::cout << "Введите данные для функции Rexlex относительно прямой\n";
		std::cout << "Выберите способ задания прямой:\n1) По двум точкам\n2) По A,B,C\n3) По k и b\n";
		int n;
		std::cin >> n;
		while (n < 1 || n > 3) {
			std::cout << "Ошибка ввода. Попробуйте ещё раз: ";
			std::cin >> n;
		}
		if (n == 1) {
			Point p1, p2;
			std::cin >> p1 >> p2;
			Line l(p1, p2);
			Selected_Shape->reflex(l);
		}
		else if (n == 2) {
			std::cout << "Введите 3 числа:\n";
			double A, B, C;
			std::cin >> A >> B >> C;
			Line l(A, B, C);
			Selected_Shape->reflex(l);
		}
		else if (n == 3) {
			std::cout << "Введите 2 числа:\n";
			double k, b;
			std::cin >> k >> b;
			Line l(k, b);
			Selected_Shape->reflex(l);
		}
		std::cout << "Изменения внесены\n\n";
		Flag_Conditions = 1;
		return;
	}
}
void Math_Methods(Buttom bt) {
	if (bt.name == "Info") {
		if (Selected_Shape->Get_Type() == "Triangle") {
			Triangle& wr = dynamic_cast<Triangle&>(*Selected_Shape);
			std::cout << "Perimetr: " << Selected_Shape->perimeter() << '\n' << "Area: " << Selected_Shape->area() << '\n' << "Ортоцентр:\n" << wr.orthocenter() << "\nЦентр масс:\n" << wr.centroid();
			
			Circle xx = wr.circumscribedCircle(), yy = wr.inscribedCircle();
			Circle* x = new Circle(xx.centre(), xx.radius()), * y = new Circle(yy.centre(), yy.radius());
			std::cout << xx << '\n' << yy;

			Shapes.push_back(x);
			Shapes.push_back(y);

			Flag_Conditions = 1;
			return;
		}
		std::cout << "Perimetr: " << Selected_Shape->perimeter() << '\n' << "Area: " << Selected_Shape->area() << std::endl << std::endl;
		Flag_Conditions = 1;
		return;
	}
	if (bt.name == "Perimetr") {
		std::cout << "Периметр выбранной фигуры: " << Selected_Shape->perimeter();
		std::cout << std::endl << std::endl;
		Flag_Conditions = 1;
		return;
	}
	if (bt.name == "  Area") {
		std::cout << "Площадь выбранной фигуры: " << Selected_Shape->area();
		std::cout << std::endl << std::endl;
		Flag_Conditions = 1;
		return;
	}
	if (bt.name == "Check_P") {
		Point p;
		std::cin >> p;
		std::string s;
		if (Selected_Shape->containsPoint(p)) s = "Yes";
		else s = "No";
		std::cout << s << std::endl << std::endl;
		Flag_Conditions = 1;
		return;
	}
	if (Flag_Conditions == 5) {
		if (flag_name == "  ==") {
			if (*Selected_Shape == *Selected_Shape_Second) {
				std::cout << "Выбранные фигуры совпадают" << std::endl << std::endl;
			}
			else std::cout << "Выбранные фигуры не совпадают" << std::endl << std::endl;
			Flag_Conditions = 1;
			return;
		}
		if (flag_name == "Similar") {
			if (Selected_Shape->isSimilarTo(*Selected_Shape_Second)) {
				std::cout << "Выбранные фигуры подобны" << std::endl << std::endl;
			}
			else std::cout << "Выбранные фигуры не подобны" << std::endl << std::endl;
			Flag_Conditions = 1;
			return;
		}
		if (flag_name == "  Cong") {
			if (Selected_Shape->isCongruentTo(*Selected_Shape_Second)) {
				std::cout << "Выбранные фигуры равны в геометрическом смысле" << std::endl << std::endl;
			}
			else std::cout << "Выбранные фигуры не равны в геометрическом смысле" << std::endl << std::endl;
			Flag_Conditions = 1;
			return;
		}
	}
	if (bt.name == "  ==" || bt.name == "Similar" || bt.name == "  Cong") {
		Flag_Conditions = 5;
		flag_name = bt.name;
		return;
	}
}

void Mouse(int button, int state, int x, int y) {
	glOrtho(-Max_Coordinate, Max_Coordinate, -Max_Coordinate, Max_Coordinate, 0., 1.);
	if ((button == GLUT_LEFT_BUTTON)) {
		if (state == GLUT_DOWN) {
			double x0 = (x - 400) / 40.0;
			double y0 = -(y - 400) / 40.0;
			Point wr(x0, y0);

			//Обработка математических операций
			if (Flag_Conditions == 4) {
				if (Buttoms_Math_Func[0].body.containsPoint(wr)) {
					Flag_Conditions = 3;
					return;
				}
				for (int i = 1; i < 8; ++i) {
					if (Buttoms_Math_Func[i].body.containsPoint(wr))
						Math_Methods(Buttoms_Math_Func[i]);
				}
			}

			if (Flag_Conditions == 5) {
				for (int i = 0; i < Shapes.size(); ++i) {
					if (Shapes[i]->containsPoint(wr)) {
						Selected_Shape_Second = Shapes[i];
						Math_Methods(Buttoms_Math_Func[0]); //В этом вызове аргумент использоваться не будет
						return;
					}
				}
			}
			

			//Обработка нажатий на фигуры
			if (Flag_Conditions == 3) {
				if (Buttoms_Func[0].body.containsPoint(wr)) {
					Flag_Conditions = 1;
					return;
				}
				if (Buttoms_Func[5].body.containsPoint(wr)) {
					Flag_Conditions = 4;
					return;
				}
				for (int i = 1; i < 7; ++i) {
					if (Buttoms_Func[i].body.containsPoint(wr))
						ReNew_Information(Buttoms_Func[i]);
				}

				return;
			}


			//Обработка кнопки New
			if (Flag_Conditions == 1) {
				Buttom NEW("New", 0);
				if (NEW.body.containsPoint(wr)) {
					Flag_Conditions = 2;
					return;
				}

				for (int i = 0; i < Shapes.size(); ++i) {
					if (Shapes[i]->containsPoint(wr)) {
						Selected_Shape = Shapes[i];
						Flag_Conditions = 3;
					}
				}
			}
			else if (Flag_Conditions == 2) {
				if (Buttoms[0].body.containsPoint(wr)) {
					Flag_Conditions = 1;
					return;
				}
				for (int i = 1; i < 7; ++i) {
					if (Buttoms[i].body.containsPoint(wr)) {
						Create_Shape(Buttoms[i].name);
						break;
					}
					Flag_Conditions = 1;
				}
			}
		}
	}
}

void Display(void) {
	glClearColor(255.0 / 255, 250.0 / 255, 205.0 / 255, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-Max_Coordinate, Max_Coordinate, -Max_Coordinate, Max_Coordinate, 0., 1.);
	

	Draw_Coordinate_Axes();
	for (int i = 0; i < Shapes.size(); ++i) {
		Shapes[i]->Draw();
	}

	switch (Flag_Conditions)
	{
	case(1):
	{
		Buttom bt("New", 0);
		bt.Draw();
		break;
	}
	case(2):
	{
		for (int i = 0; i < 7; ++i) {
			Buttoms[i].Draw();
		}
		break;
	}
	case(3):
	{
		for (int i = 0; i < 7; ++i) {
			Buttoms_Func[i].Draw();
		}
		break;
	}
	case(4):
	{
		Buttoms_Math_Func[0].R = 255;
		Buttoms_Math_Func[0].G = 0;
		Buttoms_Math_Func[0].B = 0;
		Buttoms_Math_Func[0].Draw();
		for (int i = 1; i < 8; ++i) {
			Buttoms_Math_Func[i].R = 204;
			Buttoms_Math_Func[i].G = 204;
			Buttoms_Math_Func[i].B = 255;
			Buttoms_Math_Func[i].Draw();
		}
		break;
	}
	case(5):
	{
		for (int i = 0; i < 8; ++i) {
			Buttoms_Math_Func[i].R = 204;
			Buttoms_Math_Func[i].G = 85;
			Buttoms_Math_Func[i].B = 0;
			Buttoms_Math_Func[i].Draw();
		}
		break;
	}
	default:
	{
		Flag_Conditions = 1;
		break;
	}
	}

	glutPostRedisplay();
	glutSwapBuffers();
}

int main(int argc, char** argv) {
	setlocale(LC_ALL, "Russian");
	srand(time(NULL));

	//Инициализация
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(300, 100);
	glutInitWindowSize(800, 800);
	glutCreateWindow("Window");
	glutDisplayFunc(Display);
	glutMouseFunc(Mouse);

	{
		Buttom bt1("Esc", 0, 255, 0, 0);
		Buttom bt2("Polygon", 2);
		Buttom bt3("Ellipse", 3.5);
		Buttom bt4("Circle", 5);
		Buttom bt5("Rectangle", 6.5);
		Buttom bt6("Square", 8);
		Buttom bt7("Triangle", 9.5);
		Buttoms.push_back(bt1);
		Buttoms.push_back(bt2);
		Buttoms.push_back(bt3);
		Buttoms.push_back(bt4);
		Buttoms.push_back(bt5);
		Buttoms.push_back(bt6);
		Buttoms.push_back(bt7);
	}
	{
		Buttom bt1("Esc", 0, 255, 0, 0);
		Buttom bt2("Scale", 2);
		Buttom bt3("P_Reflex", 3.5);
		Buttom bt4("L_Reflex", 5);
		Buttom bt5("Rotate", 6.5);
		Buttom bt6("Math", 8., 0, 0, 255);
		Buttom bt7("Del", 9.5, 255, 0, 0);
		Buttoms_Func.push_back(bt1);
		Buttoms_Func.push_back(bt2);
		Buttoms_Func.push_back(bt3);
		Buttoms_Func.push_back(bt4);
		Buttoms_Func.push_back(bt5);
		Buttoms_Func.push_back(bt6);
		Buttoms_Func.push_back(bt7);
	}
	{
		Buttom bt1("Esc", 0, 255, 0, 0);
		Buttom bt2("Perimetr", 2);
		Buttom bt3("  Area", 3.5);
		Buttom bt4("  ==", 5);
		Buttom bt5("  Cong", 6.5);
		Buttom bt6("Similar", 8.);
		Buttom bt7("Check_P", 9.5);
		Buttom bt8("Info", 11.);
		Buttoms_Math_Func.push_back(bt1);
		Buttoms_Math_Func.push_back(bt2);
		Buttoms_Math_Func.push_back(bt3);
		Buttoms_Math_Func.push_back(bt4);
		Buttoms_Math_Func.push_back(bt5);
		Buttoms_Math_Func.push_back(bt6);
		Buttoms_Math_Func.push_back(bt7);
		Buttoms_Math_Func.push_back(bt8);
	}
	
	Point p1(-3, -3), p2(-3, 3), p3(3, 3), p4(3, -3), pe1(-7,5), pe2(-2, 2);
	std::vector<Point> arr = { p1,p2,p3,p4 };
	Polygon* x = new Polygon(arr);
	Ellipse* y = new Ellipse(pe1, pe2, 8);
	Shapes.push_back(x);
	Shapes.push_back(y);

	glutMainLoop();
}
