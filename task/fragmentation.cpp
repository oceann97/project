#include "fragmentation.h"
#include <fstream>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_vector.h>
#include <iostream>

////Нужно переместить в с файл  иначе LNK2005
/// вектор, содержащий box-ы, являющиеся частью рабочего пространства
std::vector<Box> solution;
/// вектор, содержащий box-ы, не являющиеся частью рабочего пространства
std::vector<Box> not_solution;
/// вектор, содержащий box-ы, находящиеся на границе между "рабочим" и "нерабочим" пространством
std::vector<Box> boundary;
/// вектор, хранящий box-ы, анализируемые на следующей итерации алгоритма
std::vector<Box> temporary_boxes;
/*cilk::reducer<cilk::op_vector<Box>> solution;
cilk::reducer<cilk::op_vector<Box>> not_solution;
cilk::reducer<cilk::op_vector<Box>> boundary;
cilk::reducer<cilk::op_vector<Box>> temporary_boxes;*/


/// функции gj()
//------------------------------------------------------------------------------------------
double g1(double x1, double x2)
{
	return (x1*x1 + x2*x2 - g_l1_max*g_l1_max);
}

//------------------------------------------------------------------------------------------
double g2(double x1, double x2)
{
	return (g_l1_min*g_l1_min - x1*x1 - x2*x2);
}

//------------------------------------------------------------------------------------------
double g3(double x1, double x2)
{
	return (x1*x1 + x2*x2 - g_l2_max*g_l2_max);
}

//------------------------------------------------------------------------------------------
double g4(double x1, double x2)
{
	return (g_l2_min*g_l2_min - x1*x1 - x2*x2);
}


//------------------------------------------------------------------------------------------
low_level_fragmentation::low_level_fragmentation(double& min_x, double& min_y, double& x_width, double& y_height )
{
	current_box = Box( min_x, min_y, x_width, y_height );
}

//------------------------------------------------------------------------------------------
low_level_fragmentation::low_level_fragmentation(const Box& box)
{
	current_box = box;
}

//------------------------------------------------------------------------------------------
// Функция VerticalSplitter() разбивает переданный в качестве параметра box по ширине
void low_level_fragmentation::VerticalSplitter(const Box& box, boxes_pair& vertical_splitter_pair)
{
	double min_x, min_y, x_width, y_height;
	box.GetParameters(min_x, min_y, x_width, y_height);
	//box.GetWidhtHeight(x_width, y_height);
	double half_x_width = x_width / 2;
	Box tmp_box_1 = Box(min_x, min_y, half_x_width, y_height);
	Box tmp_box_2 = Box(min_x + half_x_width, min_y, half_x_width, y_height);
	vertical_splitter_pair = std::pair<Box, Box>(tmp_box_1, tmp_box_2);
}

//------------------------------------------------------------------------------------------
// Функция HorizontalSplitter() разбивает переданный в качестве параметра box по высоте
void low_level_fragmentation::HorizontalSplitter(const Box& box, boxes_pair& horizontal_splitter_pair)
{
	double min_x, min_y, x_width, y_height;
	box.GetParameters(min_x, min_y, x_width, y_height);
	//box.GetWidhtHeight(x_width, y_height);
	double half_y_height = y_height / 2;
	Box tmp_box_1 = Box(min_x, min_y, x_width, half_y_height);
	Box tmp_box_2 = Box(min_x, min_y + half_y_height, x_width, half_y_height);
	horizontal_splitter_pair = std::pair<Box, Box>(tmp_box_1, tmp_box_2);
}

//------------------------------------------------------------------------------------------
//Функция GetNewBoxes() разбивает переданный в качестве параметра box по большей стороне,
//вызывая VerticalSplitter() или HorizontalSplitter()
void low_level_fragmentation::GetNewBoxes(const Box& box, boxes_pair& new_pair_of_boxes)
{
	double x_width;
	double y_height;
	box.GetWidhtHeight(x_width, y_height);
	if (x_width > y_height) {
		VerticalSplitter(box, new_pair_of_boxes);
	}
	else {
		HorizontalSplitter(box, new_pair_of_boxes);
	}
}

//------------------------------------------------------------------------------------------
unsigned int low_level_fragmentation::FindTreeDepth()
{
	double box_diagonal = current_box.GetDiagonal();

	if (box_diagonal <= g_precision)
	{
		return 0;
	}
	else
	{
		boxes_pair new_boxes;
		// допустим, разобьем начальную область по ширине
		VerticalSplitter(current_box, new_boxes);
		unsigned int tree_depth = 1;

		box_diagonal = new_boxes.first.GetDiagonal();

		if (box_diagonal <= g_precision)
		{
			return tree_depth;
		}
		else
		{
			for (;;)
			{
				GetNewBoxes(new_boxes.first, new_boxes);
				++tree_depth;
				box_diagonal = new_boxes.first.GetDiagonal();

				if (box_diagonal <= g_precision)
				{
					break;
				}
			}
			return tree_depth;
		}
	}
}

//------------------------------------------------------------------------------------------
// функция ClasifyBox() анализирует box и классифицирует его
int low_level_fragmentation::ClasifyBox(const min_max_vectors& vects)
{
	//To this function we use " each box satisfies one of the following statements 4,5 or 6..."
	//"It should be noted that boxes satisfying(4) contain internal points of X; boxes satisfying(5) do not intersect with X."
	int j = 0;
	while (j < vects.first.size()) {
		double current_value = vects.first[j];
		if (current_value > 0) return 5; //Suitable for 5 condition 
		j++;
	}

	j = 0;
	while (j < vects.second.size()) {
		double current_value = vects.second[j];
		if (current_value > 0) return 6; //Suitable for 6 condition 
		j++;
	}

	return 4;//Otherwise 4 condition
}

//------------------------------------------------------------------------------------------
// Функция GetBoxType() добавляет классифицированный ранее box во множество решений, 
// или удаляет его из анализа, или добавляет его к граничной области, 
// или относит его к тем, что подлежат дальнейшему анализу
void low_level_fragmentation::GetBoxType(const Box& box)
{
	min_max_vectors vects;
	GetMinMax(box, vects);
	int box_type = ClasifyBox(vects);
	switch (box_type) {
	case 4:
		solution.push_back(box);
		break;
	case 5:
		not_solution.push_back(box);
		break;
	case 6:
	{
		boxes_pair new_boxes;
		GetNewBoxes(box, new_boxes);//Create new boxes from box
		if (new_boxes.first.GetDiagonal() < g_precision) {
			//Add to boundary
			boundary.push_back(new_boxes.first);
			boundary.push_back(new_boxes.second);
		}
		else {
			//Add to futher analysis
			temporary_boxes.push_back(new_boxes.first);
			temporary_boxes.push_back(new_boxes.second);
		}
		break;
	}
	default:
		break;
	}
}


//------------------------------------------------------------------------------------------
high_level_analysis::high_level_analysis( double& min_x, double& min_y, double& x_width, double& y_height ) :
					low_level_fragmentation(min_x, min_y, x_width, y_height) {}

//------------------------------------------------------------------------------------------
high_level_analysis::high_level_analysis( Box& box ) : low_level_fragmentation( box ) {}

//------------------------------------------------------------------------------------------
void high_level_analysis::GetMinMax( const Box& box, min_max_vectors& min_max_vecs )
{
	std::vector<double> g_min;
	std::vector<double> g_max;

	double a1min, a2min, a1max, a2max;
	double xmin, xmax, ymin, ymax;

	box.GetParameters(xmin, ymin, xmax, ymax);

	xmax = xmin + xmax;
	ymax = ymin + ymax;

	double curr_box_diagonal = box.GetDiagonal();

	if (curr_box_diagonal <= g_precision)
	{
		g_min.push_back(0);
		g_max.push_back(0);

		min_max_vecs.first = g_min;
		min_max_vecs.second = g_max;

		return;
	}

	// MIN
	// функция g1(x1,x2)
	a1min = __min(abs(xmin), abs(xmax));
	a2min = __min(abs(ymin), abs(ymax));
	g_min.push_back(g1(a1min, a2min));

	// функция g2(x1,x2)
	a1min = __max(abs(xmin), abs(xmax));
	a2min = __max(abs(ymin), abs(ymax));
	g_min.push_back(g2(a1min, a2min));

	// функция g3(x1,x2)
	a1min = __min(abs(xmin - g_l0), abs(xmax - g_l0));
	a2min = __min(abs(ymin), abs(ymax));
	g_min.push_back(g3(a1min, a2min));

	// функция g4(x1,x2)
	a1min = __max(abs(xmin - g_l0), abs(xmax - g_l0));
	a2min = __max(abs(ymin), abs(ymax));
	g_min.push_back(g4(a1min, a2min));

	// MAX
	// функция g1(x1,x2)
	a1max = __max(abs(xmin), abs(xmax));
	a2max = __max(abs(ymin), abs(ymax));
	g_max.push_back(g1(a1max, a2max));

	// функция g2(x1,x2)
	a1max = __min(abs(xmin), abs(xmax));
	a2max = __min(abs(ymin), abs(ymax));
	g_max.push_back(g2(a1max, a2max));

	// функция g3(x1,x2)
	a1max = __max(abs(xmin - g_l0), abs(xmax - g_l0));
	a2max = __max(abs(ymin), abs(ymax));
	g_max.push_back(g3(a1max, a2max));

	// функция g4(x1,x2)
	a1max = __min(abs(xmin - g_l0), abs(xmax - g_l0));
	a2max = __min(abs(ymin), abs(ymax));
	g_max.push_back(g4(a1max, a2max));

	min_max_vecs.first = g_min;
	min_max_vecs.second = g_max;
}

//------------------------------------------------------------------------------------------
//  Функция GetSolution() запускает алгоритм нахождения рабочей области манипулятора
void high_level_analysis::GetSolution()
{
	double diagonal = current_box.GetDiagonal();
	std::cout << "diagonal: "<< diagonal << std::endl;
	if (diagonal <= g_precision) solution.push_back(current_box);
	else {
		temporary_boxes.push_back(current_box);
		while (diagonal > g_precision) {
			int node_size = temporary_boxes.size();
			//std::cout << "node size: " << node_size << std::endl;
			cilk_for (int j = 0; j < node_size; j++) {
				GetBoxType(temporary_boxes.at(j));
			}
			std::cout << "temporary boxes number: "<<temporary_boxes.size() << std::endl;
			temporary_boxes.erase(temporary_boxes.begin(), temporary_boxes.begin() + node_size);
			if (!temporary_boxes.empty()) {
				diagonal = temporary_boxes.at(0).GetDiagonal();
				std::cout << "diagonal: "<< diagonal << std::endl;
			}
			else 
				return;
		}
	}
}

/*cilk::reducer<cilk::op_vector<Box>> temporary_boxes_par;

void high_level_analysis::GetSolution()
{
	double diagonal = current_box.GetDiagonal();
	std::cout << "diagonal: " << diagonal << std::endl;
	if (diagonal <= g_precision) solution.push_back(current_box);
	else {
		temporary_boxes.push_back(current_box);
		while (diagonal > g_precision) {
			int node_size = temporary_boxes.size();
			//std::cout << "node size: " << node_size << std::endl;
			temporary_boxes_par.move_in(temporary_boxes) ;
			cilk_for (int j = 0; j < node_size; j++) {
				GetBoxType(temporary_boxes_par);
			}
			temporary_boxes_par.move_out(temporary_boxes);
			std::cout << "temporary boxes number: " << temporary_boxes.size() << std::endl;
			temporary_boxes.erase(temporary_boxes.begin(), temporary_boxes.begin() + node_size);
			if (!temporary_boxes.empty()) {
				diagonal = temporary_boxes.at(0).GetDiagonal();
				std::cout << "diagonal: " << diagonal << std::endl;
			}
			else
				return;
		}
	}
}*/

//------------------------------------------------------------------------------------------
// Функция WriteResults() записывает параметры полученных box-ов (относящихся к рабочему
// пространству, к граничной области и ко множеству, не являющемуся решением) в выходные 
//файлы для дальнейшей визуализации
void WriteResults( const char* file_names[] )
{
	std::ofstream out;    // Поток записи
	out.open(file_names[0]); // Открытие файла 0 для записи
	double x_min, y_min, width, height;
	Box tmp_box;
	if (out.is_open())
	{
		for (int i = 0; i < solution.size(); i++)
		{
			tmp_box = solution.at(i); //вектор, содержащий box - ы, являющиеся частью рабочего пространства
			tmp_box.GetParameters(x_min, y_min, width, height);
			out << x_min << " " << y_min << " " << width << " " << height << std::endl;
		}
	}

	std::cout << std::endl;
	out.close();//Закрыть файловый поток

	out.open(file_names[1]); // Открытие файла 1 для записи
	if (out.is_open())
	{
		for (int i = 0; i < not_solution.size(); i++)
		{
			tmp_box = not_solution.at(i);// вектор, содержащий box-ы, не являющиеся частью рабочего пространства
			tmp_box.GetParameters(x_min, y_min, width, height);
			out << x_min << " " << y_min << " " << width << " " << height << std::endl;
		}
	}
	std::cout << std::endl;
	out.close();//Закрыть файловый поток

	out.open(file_names[2]); // Открытие файла 2 для записи
	if (out.is_open())
	{
		for (int i = 0; i < boundary.size(); i++)
		{
			tmp_box = boundary.at(i);// вектор, содержащий box-ы, находящиеся на границе между "рабочим" и "нерабочим" пространством
			tmp_box.GetParameters(x_min, y_min, width, height);
			out << x_min << " " << y_min << " " << width << " " << height << std::endl;
		}
	}
	std::cout << std::endl;
	out.close();//Закрыть файловый поток

}