#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>

using namespace std;

double b = 1.0;
double pi = acos(-1.0);
double W[3];
double eps = 0.1;
FILE *f, *f_error, *f_sigmoid, *f_anim;


struct Point {
	double x1;
	double x2;

	int target;
};
vector<Point> train_points;


double double_rand() {
  return double(rand()) / (double(RAND_MAX) + 1.0);
}


double sigmoid(double u){
	return 1.0 / (1.0 + exp(-b * u));
}


double derivative_sigmoid(double u){
	return b * sigmoid(u) * (1.0 - sigmoid(u));
}


double weighted_sum(double W[3], Point p){
	double sum = 0.0;
	sum += W[0] * 1.0;
	sum += W[1] * p.x1;
	sum += W[2] * p.x2;
	return sum;
}


double error_func(vector<Point> points){
	double error = 0.0;
	for(Point p: points){
		if (sigmoid(weighted_sum(W, p)) >= 0.5) {
			error += pow((1.0 - p.target), 2.0);
		} else {
			error += pow((0.0 - p.target), 2.0);
		}
	}
	return error *= 0.5;
}


double dE_dw0(vector<Point> points){
	double sum = 0.0;
	for (Point p: points){
		sum += (sigmoid(weighted_sum(W, p)) - p.target) * derivative_sigmoid(weighted_sum(W, p)) * 1.0;
	}
	return sum;
}


double dE_dw1(vector<Point> points){
	double sum = 0.0;
	for (Point p: points){
		sum += (sigmoid(weighted_sum(W, p)) - p.target) * derivative_sigmoid(weighted_sum(W, p)) * p.x1;
	}
	return sum;
}


double dE_dw2(vector<Point> points){
	double sum = 0.0;
	for (Point p: points){
		sum += (sigmoid(weighted_sum(W, p)) - p.target) * derivative_sigmoid(weighted_sum(W, p)) * p.x2;
	}
	return sum;
}


void make_ellipse(double x0, double y0, double angle, double a, double b, vector<Point> &points, int target, int color){
	double x, y;
	double ellipse_area;

	double angle_rad = pi / 180.0  * (angle + 90.0);

	for (x = -10.0; x < 50.0; x+= 2.5){
		for (y = -10.0; y < 50.0; y+=2.5){

			double _x = x;// + double_rand();
			double _y = y;// + double_rand();

			ellipse_area = (
				((_x - x0)*cos(angle_rad) + (_y - y0)*sin(angle_rad)) * ((_x - x0)*cos(angle_rad) + (_y - y0)*sin(angle_rad)) / a / a
				+
				((-(_x - x0)*sin(angle_rad)) + (_y - y0)*cos(angle_rad)) * ((-(_x - x0)*sin(angle_rad)) + (_y - y0)*cos(angle_rad)) / b / b
			);

			if (ellipse_area <= 1.0) {
				fprintf(f, "%lf %lf %d\n", _x, _y, color);
				points.push_back(Point {_x, _y, target});
			}
		}
	}
}


void fit(vector<Point> points){
	int iter = 0;
	int epoch = 0;
	double nu = 1.0;

	char c;

	double de_dw0 = 0.0;
	double de_dw1 = 0.0;
	double de_dw2 = 0.0;

	while (abs(error_func(points)) > eps) {
		double x,y;
		int count = 0;
		for (x = -10.0; x < 50.0; x += 0.1){
			for (y = -10.0; y < 50.0; y += 0.1){
				if (count % 19 == 0)
					;//fprintf(f_sigmoid, "%lf %lf %lf\n", x, y, sigmoid(weighted_sum(W, Point {x, y, -1})));
				count ++;
			}
		}
		for (Point p: points){
			fprintf(f_sigmoid, "%lf %lf %d\n", p.x1, p.x2, p.target+2);
		}
		fprintf(f_sigmoid, "%lf %lf 0\n", 0.0, 0.0);
		fprintf(f_sigmoid, "\n");
		for (Point p: points) {
			iter ++;

			de_dw0 = 1.0 * (sigmoid(weighted_sum(W, p)) - p.target );
			de_dw1 = p.x1 * (sigmoid(weighted_sum(W, p)) - p.target );
			de_dw2 = p.x2 * (sigmoid(weighted_sum(W, p)) - p.target );

			//cout << "de_dw" << endl;
			//cout << de_dw0 << " " << de_dw1 << " " << de_dw2 << endl;

			W[0] -= nu * de_dw0;
			W[1] -= nu * de_dw1;
			W[2] -= nu * de_dw2;

			cout << "wieghts" << endl;
			for (double w: W){
				cout << w << endl;
			}

			cout << "error local " << p.target << endl;
			cout << sigmoid(weighted_sum(W, p)) - p.target  << endl;

			cout << "Error total" << endl;
			cout << error_func(points) << endl;

			de_dw0 = 0.0;
			de_dw1 = 0.0;
			de_dw2 = 0.0;
		}

			//cin >> c;
		
		epoch ++;
		cout << "epoch " << epoch << endl;
		fprintf(f_error, "%d %lf\n", epoch, error_func(points));
	}
	fprintf(f_anim, "do for[a=1:%d]{\n\tpl 'sigmoid_3d.txt' u 1:2:3 every :::a::a w p pt 7 palette\n\tpause 0.1\n}\npause -1\n", epoch);
}


int main(int argc, char* argv[]){
	srand (time(NULL));

  	if ((f = fopen ("graph_data.txt","w")) == NULL) {
    	fprintf(stderr, "Can't open file to write\n");
  	}

  	if ((f_error = fopen ("graph_error.txt","w")) == NULL) {
    	fprintf(stderr, "Can't open file to write\n");
  	}

  	if ((f_sigmoid = fopen ("sigmoid_3d.txt","w")) == NULL) {
    	fprintf(stderr, "Can't open file to write\n");
  	}

  	if ((f_anim = fopen ("anim.gp","w")) == NULL) {
    	fprintf(stderr, "Can't open file to write\n");
  	}

	make_ellipse(25.0, 15.0, 10.0, 5.0, 15.0, train_points, 0, 2);
	//make_ellipse(25.0, 25.0, 45.0, 7.0, 15.0, train_points, 0, 2);
	fprintf(f, "\n\n");
	//make_ellipse(30.0, 30.0, 10.0, 5.0, 15.0, train_points, 1, 1);
	make_ellipse(13.0, 33.0, 70.0, 7.0, 15.0, train_points, 1, 1);
	//make_ellipse(0, 15.0, 0.0, 5.0, 5.0, train_points, 1, 1);
	fprintf(f, "\n\n");
	fprintf(f, "%lf %lf 0\n", 50.0, 50.0);
	fprintf(f, "%lf %lf 0\n", -10.0, -10.0);

	//train_points.push_back({25.0, 29.0, 1});

	/*train_points.push_back({0.0, 0.0, 0});
	train_points.push_back({0.9, 0.0, 0});
	train_points.push_back({0.0, 0.9, 0});
	train_points.push_back({0.5, 0.4, 0});
	train_points.push_back({0.4, 0.5, 0});
	train_points.push_back({1.0, 1.0, 1});
	train_points.push_back({1.1, 0.0, 1});
	train_points.push_back({0.0, 1.1, 1});
	train_points.push_back({0.6, 0.5, 1});
	train_points.push_back({0.5, 0.6, 1});*/

	cout << "points" << endl;
	for (Point p: train_points){
		cout << p.x1 << " " << p.x2 << " " << p.target << endl;
	}

	W[0] = double_rand(); // [-80 : -50]
	W[1] = double_rand(); // [-5 : 5]
	W[2] = double_rand(); //[0 : 10]

	//W[0] = 0.001;
	//W[1] = -0.5;
	//W[2] = 0.3;

	cout << "W" << endl;
	for (double w: W){
		cout << w << endl;
	}

	cout << "WEIGHT SUM" << endl;
	cout << sigmoid(weighted_sum(W, train_points[30])) << endl;

	cout << "ERROR" << endl;
	cout << error_func(train_points) << endl;
	cout << dE_dw0(train_points) << endl;
	cout << dE_dw1(train_points) << endl;
	cout << dE_dw2(train_points) << endl;

	fit(train_points);

	cout << "Point (40.0 4.0) is " << sigmoid(weighted_sum(W, Point {40.0, 4.0, -1})) << endl;
	cout << "Point (4.0 40.0) is " << sigmoid(weighted_sum(W, Point {4.0, 40.0, -1})) << endl;

	fprintf(f, "%lf %lf 3\n", 40.0, 4.0);
	fprintf(f, "%lf %lf 0\n", 4.0, 40.0);

	double x, y;

	for (x = 0.0; x <50.0; x += 0.1){
			for (y = 0.0; y < 50.0; y += 0.1){
				fprintf(f_sigmoid, "%lf %lf %lf\n", x, y, sigmoid(weighted_sum(W, Point {x, y, -1})));
			}
		}
	for (Point p: train_points){
		fprintf(f_sigmoid, "%lf %lf %d\n", p.x1, p.x2, p.target+2);
	}
	fprintf(f_sigmoid, "\n");

	FILE *w_err;

	if ((w_err = fopen ("w_error.txt","w")) == NULL) {
    	fprintf(stderr, "Can't open file to write\n");
  	}

	return 0;
}