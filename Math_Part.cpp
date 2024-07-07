//#include <iostream>
//#include <vector>
//#include <set>
//#include <cmath>
//#include <algorithm>
//#include <list>
//
//const double PI = acos(-1);
//
//struct Point {
//	//Fields
//	double x;
//	double y;
//
//	//Constructors
//	Point() = default;
//	Point(double x, double y) {
//		this->x = x;
//		this->y = y;
//	}
//	Point(std::pair<double, double> other) {
//		x = other.first;
//		y = other.second;
//	}
//
//	//Methods
//	double Dist_Between_Points(const Point& p) {
//		return sqrt((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y));
//	}
//
//	//Operators
//	void operator =(const Point& other) {
//		this->x = other.x;
//		this->y = other.y;
//	}
//	void operator =(const std::pair<double, double> other) {
//		x = other.first;
//		y = other.second;
//	}
//
//	//Bool operators
//	bool operator ==(const Point& other) const  {
//		return (this->x == other.x) && (this->y == other.y);
//	}
//	bool operator !=(const Point& other) const  {
//		return (this->x != other.x) || (this->y != other.y);
//	}
//};
//std::ostream& operator <<(std::ostream& cout, const Point& other) {
//	std::cout << "X:" << other.x << " Y:" << other.y;
//	return cout;
//}
//std::istream& operator >>(std::istream& cin, Point& other) {
//	std::cout << "Введите координаты точки.\nX:";
//	std::cin >> other.x;
//	std::cout << "Y:";
//	std::cin >> other.y;
//	return cin;
//}
//
//class Line {
//private:
//	double _A;
//	double _B;
//	double _C;
//public:
//	//Constructors
//	Line() = default;
//	Line(const double& A, const double& B, const double& C) {
//		this->_A = A;
//		this->_B = B;
//		this->_C = C;
//	}
//	//По двум точкам
//	Line(const Point& p1, const Point& p2) {
//		if (p1 == p2) throw std::overflow_error("Points must be different");
//		this->_A = p2.y - p1.y;
//		this->_B = p1.x - p2.x;
//		this->_C = p2.x * p1.y - p1.x * p2.y;
//	}
//	//По угловому коэффициенту и сдвигу (y=kx+b)
//	Line(const double& k, const double& b) {
//		this->_A = -k;
//		this->_B = 1;
//		this->_C = -b;
//	}
//	//По точке и угловому коэффициенту
//	Line(const Point& p, const double& k) {
//		this->_A = -k;
//		this->_B = 1;
//		this->_C = -(_A * p.x + p.y);
//	}
//
//	//Getters
//	double GetA() { return _A; }
//	double GetB() { return _B; }
//	double GetC() { return _C; }
//
//	//Setters
//	void Set_A(double new_A) { _A = new_A; }
//	void Set_B(double new_B) { _B = new_B; }
//	void Set_C(double new_C) { _C = new_C; }
//
//	//Methods
//	Point  point_of_intersection(Line& other) {
//		Point ans;
//		ans.x = (other.GetB() - _B + other.GetC() - _C) / (_A - other.GetA());
//		ans.y = (other.GetA() - _A + other.GetC() - _C) / (_B - other.GetB());
//		return ans;
//	}
//
//	//Operators
//	void operator =(const Line& other) {
//		this->_A = other._A;
//		this->_B = other._B;
//		this->_C = other._C;
//	}
//
//	//Bool operators
//	bool operator ==(const Line& other) const {
//		return (this->_A == other._A) && (this->_B == other._B) && (this->_C == other._C);
//	}
//	bool operator !=(const Line& other) const {
//		return !(*this == other);
//	}
//};
//std::istream& operator >>(std::istream& cin, Line& other) {
//	double A, B, C;
//	std::cout << "Введите коэффициенты из общего уравнения прямой (Ax+By+C=0).\nA:";
//	std::cin >> A;
//	std::cout << "B:";
//	std::cin >> B;
//	std::cout << "C:";
//	std::cin >> C;
//	other.Set_A(A); other.Set_B(B); other.Set_C(C);
//	return cin;
//}
//std::ostream& operator <<(std::ostream& cout, Line& other) {
//	std::cout << "A:" << other.GetA() << " B:" << other.GetB() << " C:" << other.GetC();
//	return cout;
//}
//
////-----------------------------------------------------------------------------------------------------------------------//
//
//class Shape {
//public:
//	//Fields
//	std::string Type = "Shape";
//
//	//Getters
//	std::string Get_Type() {
//		return Type;
//	}
//
//	//glut
//	virtual void Draw() = 0;
//	virtual void rotate(Point center, double angle) = 0;
//	virtual void reflex(Line axis) = 0;
//	virtual void reflex(Point center) = 0;
//	virtual void scale(Point center, double coef) = 0;
//
//	//Math-functions
//	virtual double  perimeter() = 0;
//	virtual bool operator ==(Shape& other) = 0;
//	virtual bool isCongruentTo(Shape& other) = 0;
//	virtual double area() = 0;
//	virtual bool containsPoint(Point& other) = 0;
//	virtual bool isSimilarTo(Shape& other) = 0;
//};
//
//class Polygon : public Shape {
//protected:
//	std::vector<Point> _vertices;
//	int _verticesCount;
//	std::string Type = "Polygon";
//public:
//	//Constructors
//	Polygon() = default;
//	Polygon(std::vector<Point> vertices) {
//		_verticesCount = vertices.size();
//		_vertices.resize(_verticesCount);
//		for (int i = 0; i < _verticesCount; ++i) {
//			_vertices[i] = vertices[i];
//		}
//	}
//	Polygon(std::initializer_list<Point> list) {
//		for (auto& x : list) {
//			_vertices.push_back(x);
//		}
//		_verticesCount = _vertices.size();
//	}
//	
//	//Getters
//	int Get_verticesCount() { 
//		return _verticesCount; 
//	}
//	const std::vector<Point>& Get_vertices() {
//		return _vertices;
//	}
//
//	//Methods
//	bool isConvex() {
//		if (_verticesCount == 3) return true;
//		for (int i = 0; i < _verticesCount - 1; ++i) {
//			Line wr_line(_vertices[i], _vertices[i + 1]);
//			for (int j = i + 2; j < _verticesCount - 1; ++j) {
//				double fx1; // Variable to store a * x1 + b * y1 - c
//				double fx2; // Variable to store a * x2 + b * y2 - c
//
//				fx1 = wr_line.GetA() * _vertices[j].x + wr_line.GetB() * _vertices[j].y - wr_line.GetC();
//				fx2 = wr_line.GetA() * _vertices[j + 1].x + wr_line.GetB() * _vertices[j + 1].y - wr_line.GetC();
//
//				if (fx1 * fx2 < 0) return false;
//			}
//		}
//		return true;
//	}
//
//	//Override glut-methods
//	void Draw() override {}
//	void rotate(Point center, double angle) override {
//		std::vector<Point> v(_vertices.size());
//		for (int i = 0; i < _vertices.size(); ++i) {
//			v[i].x = center.x + (_vertices[i].x - center.x) * cos(angle / 180 * PI) - (_vertices[i].y - center.y) * sin(angle / 180 * PI);
//			v[i].y = center.y + (_vertices[i].x - center.x) * sin(angle / 180 * PI) + (_vertices[i].y - center.y) * cos(angle / 180 * PI);
//		}
//		_vertices = v;
//	}
//	void reflex(Line axis) override {
//		int sz = _vertices.size();
//		std::vector<Point> v = _vertices;
//		if (axis.GetA() == 0) {
//			for (int i = 0; i < sz; ++i) {
//				_vertices[i].y = 2 * (-axis.GetC() / axis.GetB()) - _vertices[i].y;
//			}
//		}
//		else if (axis.GetB() == 0) {
//			for (int i = 0; i < sz; ++i) {
//				_vertices[i].x = 2 * (-axis.GetC() / axis.GetB()) - _vertices[i].x;
//			}
//		}
//		else {
//			double k = -axis.GetA() / axis.GetB(), c = -axis.GetC() / axis.GetB();
//			for (int i = 0; i < sz; ++i) {
//				v[i].x = 2 * (_vertices[i].x - k * (c - _vertices[i].y)) / (1 + k * k) - _vertices[i].x;
//				v[i].y = (k * (2 * _vertices[i].x + k * _vertices[i].y) + 2 * c - _vertices[i].y) / (1 + k * k);
//			}
//			_vertices = v;
//		}
//	}
//	void reflex(Point center) override {
//		for (int i = 0; i < _vertices.size(); ++i) {
//			_vertices[i].x += 2 * (center.x - _vertices[i].x);
//			_vertices[i].y += 2 * (center.y - _vertices[i].y);
//		}
//	}
//	void scale(Point center, double coef) override {
//		for (int i = 0; i < _vertices.size(); ++i) {
//			_vertices[i].x = coef * (_vertices[i].x - center.x) + center.x;
//			_vertices[i].y = coef * (_vertices[i].y - center.y) + center.y;
//		}
//	}
//
//	//Override math-methods
//	double perimeter() override {
//		double ans = 0;
//		for (int i = 0; i < _verticesCount - 1; ++i) {
//			ans += _vertices[i].Dist_Between_Points(_vertices[i + 1]);
//		}
//		ans += _vertices[0].Dist_Between_Points(_vertices[_verticesCount - 1]);
//		return ans;
//	}
//	bool operator ==(Shape& other) override {
//		if (other.Type == "Ellipse" || other.Type == "Circle") {
//			return false;
//		}
//		Polygon& wr = dynamic_cast<Polygon&>(other);
//		if (_verticesCount != wr._verticesCount) return false;
//		for (int i = 0; i < _verticesCount; ++i) {
//			if (_vertices[i] != wr._vertices[i]) return false;
//		}
//		return true;
//	}
//	bool isCongruentTo(Shape& other) override {
//		if (other.Type == "Ellipse" || other.Type == "Circle") {
//			return false;
//		}
//		Polygon& wr = dynamic_cast<Polygon&>(other);
//		if (_verticesCount != wr._verticesCount) return false;
//
//		std::vector<double> arr1, arr2;
//		for (int i = 0; i < _verticesCount - 1; ++i) {
//			arr1.push_back(_vertices[i].Dist_Between_Points(_vertices[i + 1]));
//			arr2.push_back(wr._vertices[i].Dist_Between_Points(wr._vertices[i + 1]));
//		}
//		arr1.push_back(_vertices[0].Dist_Between_Points(_vertices[_verticesCount - 1]));
//		arr2.push_back(wr._vertices[0].Dist_Between_Points(wr._vertices[wr._verticesCount - 1]));
//		std::sort(arr1.begin(), arr1.end());
//		std::sort(arr2.begin(), arr2.end());
//		for (int i = 0; i < _verticesCount; ++i) {
//			if (arr1[i] != arr2[i]) return false;
//		}
//		return true;
//	}
//	double area() override {
//		int sz = _verticesCount;
//		double s = 0;
//		for (int i = 0; i < sz - 1; ++i) {
//			s += _vertices[i].x * _vertices[i + 1].y;
//		}
//		s += _vertices[sz - 1].x * _vertices[0].y;
//		for (int i = 0; i < sz - 1; ++i) {
//			s -= _vertices[i + 1].x * _vertices[i].y;
//		}
//		s -= _vertices[0].x * _vertices[sz - 1].y;
//		s = abs(s);
//		s /= 2;
//		return s;
//	}
//	bool containsPoint(Point& other) override {
//		//Если луч, выпущенный из точки пересекает четное количество сторон, то
//		//точка не лежит внутри многоугольника. Иначе точка лежит в многоугольнике
//		bool res = false;
//		int j = _verticesCount - 1;
//		for (int i = 0; i < _vertices.size(); i++) {
//			if ((_vertices[i].y < other.y && _vertices[j].y >= other.y || _vertices[j].y < other.y && _vertices[i].y >= other.y) &&
//				(_vertices[i].x + (other.y - _vertices[i].y) / (_vertices[j].y - _vertices[i].y) * (_vertices[j].x - _vertices[i].x) < other.x))
//				res = !res;
//			j = i;
//		}
//		return res;
//	}
//	bool isSimilarTo(Shape& other) override {
//		if (other.Type == "Ellipse" || other.Type == "Circle") return false;
//
//		Polygon& a = dynamic_cast<Polygon&>(other);
//
//		if (_vertices.size() != a._vertices.size()) return false;
//		if (_vertices.size() == 1) return true;
//		int sz = _verticesCount;
//
//		double x1 = _vertices[0].x - _vertices[sz - 1].x;
//		double y1 = _vertices[0].y - _vertices[sz - 1].y;
//		double x2 = _vertices[1].x - _vertices[0].x;
//		double y2 = _vertices[1].y - _vertices[0].y;
//		std::vector<double> vector_of_skalyars(sz);
//		std::vector<double> vector_of_lengths(sz);
//		vector_of_skalyars[0] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
//		vector_of_lengths[0] = sqrt(x1 * x1 + y1 * y1);
//		for (int i = 2; i < sz; i++) {
//			x1 = x2;
//			y1 = y2;
//			x2 = _vertices[i].x - _vertices[i - 1].x;
//			y2 = _vertices[i].y - _vertices[i - 1].y;
//			vector_of_skalyars[i - 1] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
//			vector_of_lengths[i - 1] = sqrt(x1 * x1 + y1 * y1);
//		}
//
//		x1 = x2;
//		y1 = y2;
//		x2 = _vertices[0].x - _vertices[sz - 1].x;
//		y2 = _vertices[0].y - _vertices[sz - 1].y;
//		vector_of_skalyars[sz - 1] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
//		vector_of_lengths[sz - 1] = sqrt(x1 * x1 + y1 * y1);
//
//		x1 = a._vertices[0].x - a._vertices[sz - 1].x;
//		y1 = a._vertices[0].y - a._vertices[sz - 1].y;
//		x2 = a._vertices[1].x - a._vertices[0].x;
//		y2 = a._vertices[1].y - a._vertices[0].y;
//		std::vector<double> vector_of_skalyars_sec(sz);
//		std::vector<double> vector_of_lengths_sec(sz);
//		vector_of_skalyars_sec[0] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
//		vector_of_lengths_sec[0] = sqrt(x1 * x1 + y1 * y1);
//		for (int i = 2; i < sz; i++) {
//			x1 = x2;
//			y1 = y2;
//			x2 = a._vertices[i].x - a._vertices[i - 1].x;
//			y2 = a._vertices[i].y - a._vertices[i - 1].y;
//			vector_of_skalyars_sec[i - 1] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
//			vector_of_lengths_sec[i - 1] = sqrt(x1 * x1 + y1 * y1);
//		}
//
//		x1 = x2;
//		y1 = y2;
//		x2 = a._vertices[0].x - a._vertices[sz - 1].x;
//		y2 = a._vertices[0].y - a._vertices[sz - 1].y;
//		vector_of_skalyars_sec[sz - 1] = (x1 * x2 + y1 * y2) / sqrt(x1 * x1 + y1 * y1) / sqrt(x2 * x2 + y2 * y2);
//		vector_of_lengths_sec[sz - 1] = sqrt(x1 * x1 + y1 * y1);
//
//
//		bool flag = true, id;
//		double k;
//		for (int i = 0; i < sz; i++) {
//			flag = false;
//			if (vector_of_skalyars_sec[0] == vector_of_skalyars[i]) {
//				flag = true;
//				id = i;
//				k = vector_of_lengths_sec[0] / vector_of_lengths[i];
//				for (int j = i + 1, cnt = 1; cnt < sz; ++cnt, ++j) {
//					if (j == sz) j = 0;
//					if (vector_of_skalyars_sec[cnt] != vector_of_skalyars[j] || k != (vector_of_lengths_sec[cnt] / vector_of_lengths[j])) {
//						flag = false;
//						break;
//					}
//				}
//			}
//			if (flag) break;
//		}
//
//		return flag;
//	}
//};
//std::ostream& operator <<(std::ostream& cout, Polygon& other) {
//	std::cout << "Координаты вершин многоугольника в порядке обхода:\n";
//	for (int i = 0; i < other.Get_verticesCount(); ++i) {
//		std::cout << i + 1 << ") " << other.Get_vertices()[i] << '\n';
//	}
//	return cout;
//}
//
//class Ellipse : public Shape {
//protected:
//	std::pair<Point, Point> _focuses;
//	double _a, _b, _c;
//	std::string Type = "Ellipse";
//public:
//	//Constructors
//	Ellipse() = default;
//	Ellipse(const Point& f1, const Point& f2, const double& dist) {
//		_focuses.first = f1;
//		_focuses.second = f2;
//		_a = dist / 2;
//		_c = sqrt(pow(_focuses.first.x - _focuses.second.x, 2) + pow(_focuses.first.y - _focuses.second.y, 2)) / 2;
//		_b = sqrt(_a * _a - _c * _c);
//	}
//
//	//Getters
//	std::pair<Point, Point> Get_focuses() {
//		return _focuses;
//	}
//	double Get_sumOfDistanse() {
//		return _a * 2.;
//	}
//	std::string Get_Type() {
//		return Type;
//	}
//
//	//Operators
//	void operator =(const Ellipse& other) {
//		_focuses = other._focuses;
//		_a = other._a;
//		_b = other._b;
//		_c = other._c;
//	}
//
//	//Methods
//	Point centre() {
//		Point ans;
//		ans.x = (_focuses.first.x + _focuses.second.x) / 2;
//		ans.y = (_focuses.first.y + _focuses.second.y) / 2;
//		return ans;
//	}
//	double eccentrisity() {
//		return _c / _a;
//	}
//	std::pair<Line, Line> directrices() {
//		return std::make_pair(Line(0, _b / (_c / _a)), Line(0, -_b / (_c / _a)));
//	}
//
//	//Override glut-methods
//	void Draw() override {}
//	void rotate(Point center, double angle) override {
//		Point g1, g2;
//		Point f1 = _focuses.first, f2 = _focuses.second;
//		g1.x = center.x + (f1.x - center.x) * cos(angle / 180 * PI) - (f1.y - center.y) * sin(angle / 180 * PI);
//		g1.y = center.y + (f1.x - center.x) * sin(angle / 180 * PI) + (f1.y - center.y) * cos(angle / 180 * PI);
//
//		g2.x = center.x + (f2.x - center.x) * cos(angle / 180 * PI) - (f2.y - center.y) * sin(angle / 180 * PI);
//		g2.y = center.y + (f2.x - center.x) * sin(angle / 180 * PI) + (f2.y - center.y) * cos(angle / 180 * PI);
//
//		f1 = g1;
//		f2 = g2;
//	}
//	void reflex(Point center) override {
//		Point f1 = _focuses.first, f2 = _focuses.second;
//		
//		f1.x += 2 * (center.x - f1.x);
//		f1.y += 2 * (center.y - f1.y);
//
//		f2.x += 2 * (center.x - f2.x);
//		f2.y += 2 * (center.y - f2.y);
//	}
//	void reflex(Line axis) override {
//		Point f1 = _focuses.first, f2 = _focuses.second;
//		if (axis.GetA() == 0) {
//			f1.y = 2 * (-axis.GetC() / axis.GetB()) - f1.y;
//			f2.y = 2 * (-axis.GetC() / axis.GetB()) - f2.y;
//		}
//		else if (axis.GetB() == 0) {
//			f1.x = 2 * (-axis.GetC() / axis.GetA()) - f1.x;
//			f2.x = 2 * (-axis.GetC() / axis.GetA()) - f2.x;
//		}
//		else {
//			double k = -axis.GetA() / axis.GetB(), c = -axis.GetC() / axis.GetB();
//			Point g1, g2;
//			g1.x = 2 * (f1.x - k * (c - f1.y)) / (1 + k * k) - f1.x;
//			g1.y = (k * (2 * f1.x + k * f1.y) + 2 * c - f1.y) / (1 + k * k);
//			g2.x = 2 * (f2.x - k * (c - f2.y)) / (1 + k * k) - f2.x;
//			g2.y = (k * (2 * f2.x + k * f2.y) + 2 * c - f2.y) / (1 + k * k);
//			f1 = g1;
//			f2 = g2;
//		}
//	}
//	void scale(Point center, double coef) override {
//		Point f1 = _focuses.first, f2 = _focuses.second;
//		f1.x = coef * (f1.x - center.x) + center.x;
//		f1.y = coef * (f1.y - center.y) + center.y;
//		f2.x = coef * (f2.x - center.x) + center.x;
//		f2.y = coef * (f2.y - center.y) + center.y;
//		_a *= coef;
//		_b *= coef;
//		_c *= coef;
//	}
//
//	//Override math-methods
//	double perimeter() override {
//		return 2 * PI * sqrt((_a * _a + _b * _b) / 2);
//	}
//	bool operator ==(Shape& other) override {
//		if (other.Type != "Ellipse" && other.Type != "Circle") return false;
//		Ellipse& wr = dynamic_cast<Ellipse&>(other);
//		if (_focuses != wr._focuses || _a != wr._a || _b != wr._b || _c != wr._c) return false;
//		return true;
//	}
//	bool isCongruentTo(Shape& other) override {
//		if (other.Type != "Ellipse" && other.Type != "Circle") return false;
//		Ellipse& wr = dynamic_cast<Ellipse&>(other);
//		if (this->perimeter() != wr.perimeter()) return false;
//		if (_a != wr._a || _b != wr._b || _c != wr._c) return false;
//		return true;
//	}
//	double area() override {
//		return PI * _a * _b;
//	}
//	bool containsPoint(Point& other) override {
//		double r1, r2;
//		r1 = sqrt(pow(_focuses.first.x - other.x, 2) + pow(_focuses.first.y - other.y, 2));
//		r2 = sqrt(pow(_focuses.second.x - other.x, 2) + pow(_focuses.second.y - other.y, 2));
//		return (r1 + r2 <= 2 * _a);
//	}
//	bool isSimilarTo(Shape& other) override {
//		if (other.Type != "Ellipse" && other.Type != "Circle") return false;
//		Ellipse& an = dynamic_cast<Ellipse&>(other);
//		return (_a / an._a == _b / an._b && _a / an._a == _c / an._c);
//	}
//};
//std::ostream& operator <<(std::ostream& cout, Ellipse& other) {
//	std::cout << "Фокусы элипса и его сумма расстояний от фокусов до контура:\nF1) " << other.Get_focuses().first << "\nF2) " << other.Get_focuses().second << "\nDistanse) " << other.Get_sumOfDistanse();
//	return cout;
//}
//
//class Circle : public Ellipse {
//public:
//	//Constructors
//	Circle() : Ellipse() {
//		Type = "Circle";
//	}
//	Circle(const Point& centre, double r) : Ellipse(centre, centre, r * 2.) {
//		Type = "Circle";
//	}
//
//	//Operators
//	void operator =(const Circle& other) {
//		this->_a = other._a;
//		this->_focuses = other._focuses;
//		this->_b = other._b;
//		this->_c = other._c;
//	}
//
//	//Getters
//	Point centre() {
//		return _focuses.first;
//	}
//	double radius() {
//		return _a;
//	}
//	std::string Get_Type() {
//		return Type;
//	}
//};
//std::ostream& operator <<(std::ostream& cout, Circle& other) {
//	std::cout << "Координаты центра окружности и её радиус:\nЦентр) " << other.Get_focuses().first << "\nРадиус) " << other.radius();
//	return cout;
//}
//
//class Rectangle : public Polygon {
//public:
//	//Constructors
//	Rectangle() = default;
//	Rectangle(std::vector<Point> vertices) : Polygon(vertices) {
//		Type = "Rectangle";
//	}
//	Rectangle(std::initializer_list<Point> lst) : Polygon(lst) {
//		Type = "Rectangle";
//	}
//	Rectangle(Point p1, Point p2, double coef) {
//		Type = "Rectangle";
//
//		if (coef < 1) {
//			coef = 1 / coef;
//		}
//		Point c((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
//		Polygon pn({ p1, p2 });
//		pn.rotate(c, 180 / (coef + 1) * coef);
//		std::vector<Point> ver = pn.Get_vertices();
//		_vertices.resize(4);
//		_vertices[0] = p1;
//		_vertices[1] = ver[0];
//		_vertices[2] = p2;
//		_vertices[3] = ver[1];
//		if (_vertices[0].Dist_Between_Points(_vertices[3]) > _vertices[3].Dist_Between_Points(_vertices[2])) {
//			this->reflex(Line(p1, p2));
//		}
//		_verticesCount = 4;
//	}
//
//	//Getters
//	std::string Get_Type() {
//		return Type;
//	}
//
//	//Methods
//	std::pair<Line, Line> diagonals() {
//		Line l1(_vertices[0], _vertices[2]);
//		Line l2(_vertices[1], _vertices[3]);
//		std::pair<Line, Line> ans = { l1, l2 };
//		return ans;
//	}
//	Point centre() {
//		Rectangle copy = *this;
//		std::pair<Line, Line> diagonals = copy.diagonals();
//		return diagonals.first.point_of_intersection(diagonals.second);
//	}
//	double area() {
//		double ans = _vertices[0].Dist_Between_Points(_vertices[1]) * _vertices[1].Dist_Between_Points(_vertices[2]);
//		return ans;
//	}
//};
//
//class Square : public Rectangle {
//public:
//	//Constructors
//	Square() : Rectangle() {
//		Type = "Square";
//	}
//	Square(std::vector<Point> arr) : Rectangle(arr) {
//		Type = "Square";
//	}
//	Square(const Point& p1, const Point& p2) : Rectangle(p1, p2, 1) {
//		Type = "Square";
//	}
//	Square(std::initializer_list<Point> lst) : Rectangle(lst) {
//		Type = "Square";
//	}
//
//	//Getters
//	std::string Get_Type() {
//		return Type;
//	}
//
//	//Methods
//	Circle inscribedCircle() {
//		double wr_r = _vertices[0].Dist_Between_Points(_vertices[1]) / 2.;
//		Point wr_p((_vertices[0].x + _vertices[2].x) / 2, (_vertices[0].y + _vertices[2].y) / 2);
//		Circle ans(wr_p, wr_r);
//		return ans;
//	}
//	Circle circumscribedCircle() {
//		Point wr_p((_vertices[0].x + _vertices[2].x) / 2, (_vertices[0].y + _vertices[2].y) / 2);
//		double wr_r = _vertices[0].Dist_Between_Points(_vertices[2]) / 2;
//		Circle ans(wr_p, wr_r);
//		return ans;
//	}
//};
//
//class Triangle : public Polygon {
//	public:
//	//Constructors
//	Triangle() : Polygon() {
//		Type = "Triangle";
//	}
//	Triangle(std::vector<Point> arr) : Polygon(arr){
//		Type = "Triangle";
//	}
//	Triangle(std::initializer_list<Point> lst) : Polygon(lst) {
//		Type = "Triangle";
//	}
//
//	//Getters
//	std::string Get_Type() {
//		return Type;
//	}
//	
//	//Methods
//	Circle circumscribedCircle() {
//		Point A = _vertices[0], B = _vertices[1], C = _vertices[2];
//		double D = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
//
//		double centerX_of_Circle = ((pow(A.x, 2) + pow(A.y, 2)) * (B.y - C.y) +
//			(pow(B.x, 2) + pow(B.y, 2)) * (C.y - A.y) +
//			(pow(C.x, 2) + pow(C.y, 2)) * (A.y - B.y)) / D;
//
//		double centerY_of_Circle = ((pow(A.x, 2) + pow(A.y, 2)) * (C.x - B.x) +
//			(pow(B.x, 2) + pow(B.y, 2)) * (A.x - C.x) +
//			(pow(C.x, 2) + pow(C.y, 2)) * (B.x - A.x)) / D;
//
//		Point p(centerX_of_Circle, centerY_of_Circle);
//		double r = A.Dist_Between_Points(p);
//		return Circle(p, r);
//	}
//	Circle inscribedCircle() {
//		Point A = _vertices[0], B = _vertices[1], C = _vertices[2];
//		double a = B.Dist_Between_Points(C);
//		double b = A.Dist_Between_Points(C);
//		double c = A.Dist_Between_Points(B);
//		double wr = a + b + c;
//		
//		double centerX_of_Circle = (a * A.x + b * B.x + c * C.x) / wr;
//		double centerY_of_Circle = (a * A.y + b * B.y + c * C.y) / wr;
//		Point p(centerX_of_Circle, centerY_of_Circle);
//		
//		double s = this->area();
//		double r = s / (a + b + c / 2);
//		return Circle(p, r);
//	}
//	Point centroid() {//Точка пересечения медиан
//		Point wr_point((_vertices[0].x + _vertices[1].x) / 2, (_vertices[0].y + _vertices[1].y) / 2);
//		return Point((_vertices[2].x + 2 * wr_point.x) / 3, (_vertices[2].y + 2 * wr_point.y) / 3);
//	}
//	Point orthocenter() {//Точка пересечения высот
//		Point ans;
//		Point A = _vertices[0], B = _vertices[1], C = _vertices[2];
//		double a = B.Dist_Between_Points(C);
//		double b = A.Dist_Between_Points(C);
//		double c = A.Dist_Between_Points(B);
//		ans.x = (A.x * a * a + B.x * b * b + C.x * c * c) / (a * a + b * b + c * c);
//		ans.y = (A.y * a * a + B.y * b * b + C.y * c * c) / (a * a + b * b + c * c);
//		return ans;
//	}
//	Line EulerLine() {
//		return Line(this->centroid(), this->orthocenter());
//	}
//	Circle ninePointsCircle() {
//		Point wr1 = (this->inscribedCircle()).centre();
//		Point wr2 = this->orthocenter();
//		Point wr3((wr1.x + wr2.x) / 2, (wr1.y + wr2.y) / 2);
//		Point wr4((_vertices[0].x + _vertices[1].x) / 2, (_vertices[0].y + _vertices[1].y) / 2);
//		double r = wr3.Dist_Between_Points(wr4);
//		return Circle(wr3, r);
//	}
//};
//
////-----------------------------------------------------------------------------------------------------------------------//
//
//
////int main() {
////	setlocale(LC_ALL, "Russian");
////	
////	Point p1(-1, 1);
////	Point p2(1, 1);
////	//bool wr = (p1 == p2);
////	Point p3(1, -1);
////	Point p4(0, -1.5);
////	Point p5(-1, -1);
////	std::vector<Point> arr = {p1, p2, p3, p4, p5};
////	//Square pl1(arr);
////	Polygon pl2(arr);
////	//std::cout << pl2.area();
////	//std::cout << pl1.area();
////	Ellipse el(p1, p2, 7.2);
////	//std::cout << el.perimeter();
////	//Square rc(p1, p3);
////	//Circle cr = rc.circumscribedCircle();
////	
////	/*Circle cr(p1, 3);
////	std::cout << cr.Get_sumOfDistanse();*/
////
////	/*arr.push_back(p1);
////	arr.push_back(p2);
////	arr.push_back(p3);
////	arr.push_back(p4);
////	Polygon sq(arr);
////	std::cout << sq.isConvex();*/
////	/*Polygon pl = {p1, p2, p3, p4, p5};
////	std::cout << pl;*/
////}