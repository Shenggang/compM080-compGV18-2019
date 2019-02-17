#include <igl/readPLY.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/writePLY.h>
#include <igl/file_exists.h>
#include <imgui/imgui.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "mytools.h"


class MyContext
{
public:

	MyContext() :point_size(8), line_width(0.5), mode(0), prob(0.5), rotation(45), trans(-0.0520211,-0.000383981,-0.0109223)
	{
	
		//initial vertices and faces
		if (!igl::file_exists("../data/bun000.ply"))
		{
			std::cout << "[error] cannot locate model file at '../data/bun000.ply' \nPress any key to exit\n";
			char c;
			std::cin >> c;
			exit(1);
		}
		if (!igl::file_exists("../data/bun045.ply"))
		{
			std::cout << "[error] cannot locate model file at '../data/bun045.ply' \nPress any key to exit\n";
			char c;
			std::cin >> c;
			exit(1);
		}
		if (!igl::file_exists("../data/bun090.ply"))
		{
			std::cout << "[error] cannot locate model file at '../data/bun090.ply' \nPress any key to exit\n";
			char c;
			std::cin >> c;
			exit(1);
		}
		if (!igl::file_exists("../data/bun180.ply"))
		{
			std::cout << "[error] cannot locate model file at '../data/bun180.ply' \nPress any key to exit\n";
			char c;
			std::cin >> c;
			exit(1);
		}
		if (!igl::file_exists("../data/bun270.ply"))
		{
			std::cout << "[error] cannot locate model file at '../data/bun270.ply' \nPress any key to exit\n";
			char c;
			std::cin >> c;
			exit(1);
		}

		igl::readPLY("../data/bun000.ply", V_set_1, F_set_dummy);
		calculate_vertex_normal_flann(V_set_1, N_set_1); // Set 1 is fixed, normal need to be calculated only once
	}		
	~MyContext() {}	

	Eigen::MatrixXd V_set_1;
	Eigen::MatrixXd N_set_1;
	Eigen::MatrixXd V_set_2;
	Eigen::MatrixXd N_set_2;
	std::vector<Eigen::MatrixXd> V_set_vec;
	std::vector<Eigen::MatrixXd> N_set_vec;
	Eigen::MatrixXi F_set_dummy;

	Eigen::RowVector3d mycolors[5] = {
		Eigen::RowVector3d(1,0,0), Eigen::RowVector3d(0,1,0), Eigen::RowVector3d(0,0,1),
		Eigen::RowVector3d(0,1,1), Eigen::RowVector3d(1,0,1)
	};

	float point_size;
	float line_width;
	int mode;
	float prob;
	int noise_std;
	int rotation = 0;
	bool init_trans = false;
	bool use_normal = false;
	bool shade = false;
	Eigen::RowVector3d trans;

	std::vector<double> loss;
	std::vector<std::vector<double>> que5_loss;

	void load_five_meshes()
	{
		std::string files[5] = {"../data/bun000.ply", "../data/bun045.ply", "../data/bun090.ply", 
						"../data/bun270.ply", "../data/bun315.ply"};
		V_set_vec.clear();
		for (int i = 0; i < 5; ++i)
		{
			Eigen::MatrixXd V,N;
			igl::readPLY(files[i], V, F_set_dummy);
			V_set_vec.push_back(V);
			if (N_set_vec.size() < 5)
			{
				if (i == 0){
				calculate_vertex_normal_flann(V,N);
				}
				N_set_vec.push_back(N);
			}
		}
	}

	void initialise_point_cloud()
	{
		if (mode == 2)
		{
			load_five_meshes();
			loss.clear();
			que5_loss.clear();
			for (int i = 0; i < 4; i++)
			{
				que5_loss.push_back(std::vector<double>());
			}
			return;
		}
		Eigen::Matrix3d rotmat;
		get_rotation_Y(rotation, rotmat);
		if (mode == 0)
		{
			igl::readPLY("../data/bun045.ply", V_set_2, F_set_dummy);

			V_set_2 = (rotmat*(V_set_2.transpose())).transpose();
			if (init_trans)
			{
				V_set_2 += trans.replicate(V_set_2.rows(),1);
			}
		} else if (mode == 1)
		{
			igl::readPLY("../data/bun000.ply", V_set_2, F_set_dummy);
			V_set_2 = (rotmat*(V_set_2.transpose())).transpose();
		}
		std::default_random_engine generator;
  		std::normal_distribution<double> distribution(0, noise_std*0.0001);

		for (int i = 0; i < V_set_2.rows(); ++i)
		{
			for (int j = 0; j < V_set_2.cols(); ++j)
			{
				V_set_2(i,j) += distribution(generator);
			}
		}
		loss.clear();
	}

	Eigen::MatrixXd chooseColor(Eigen::MatrixXd color, Eigen::MatrixXd shade_color)
	{
		return shade ? shade_color : color;
	} 

	void reset_display(igl::opengl::glfw::Viewer& viewer)
	{
		viewer.data().clear();

		if (mode == 0 || mode == 1)
		{
			if (shade)
			{
				calculate_vertex_normal_flann(V_set_2, N_set_2);
			}
			// Load in the first mesh
			viewer.data().add_points(V_set_1, chooseColor(Eigen::RowVector3d(1,0,0), N_set_1));

			// Load in the second mesh and focus
			viewer.data().add_points(V_set_2, chooseColor(Eigen::RowVector3d(0,0,1), N_set_2));
			viewer.core.align_camera_center(V_set_1);
		} else
		{
			for (int i = 0; i < 5; ++i)
			{
				if (shade && (i > 0))
				{
					std::cout << "Calculating normal for mesh " << i << "/4" << std::endl;
					calculate_vertex_normal_flann(V_set_vec[i], N_set_vec[i]);
				}
				viewer.data().add_points(V_set_vec[i], chooseColor(mycolors[i], N_set_vec[i]));
			}
		}

		viewer.data().line_width = line_width;
		viewer.data().point_size = point_size;
	}

	void init_display(igl::opengl::glfw::Viewer& viewer)
	{
		initialise_point_cloud();
		reset_display(viewer);
	}

	double icp_step(Eigen::MatrixXd const & V_1, Eigen::MatrixXd & V_2)
	{
		// Actual icp step, without normal

		// subsample
		Eigen::MatrixXd V_subset;
		if (prob < 1)
		{
			std::vector<int> indices;
			subsample(V_1.rows(), prob, indices);
			V_subset.resize(indices.size(), 3);

			for (int i = 0; i < indices.size(); ++i)
			{
				V_subset.row(i) = V_1.row(indices[i]);
			}
		} else
		{
			V_subset = V_1;
		}
		
		// setup coorespondance
		Eigen::MatrixXd V_ref;
		Eigen::MatrixXd V_NP;
		register_closest_point(V_subset, V_2, V_ref, V_NP);
		
		// calculate transformation
		Eigen::Matrix3d R;
		Eigen::Vector3d t;
		double l = estimate_rotation_translation(V_ref, V_NP, R, t);
		if (mode != 2)
		{
			loss.push_back(l);
		}
		V_2 = (R*(V_2.transpose()) + t.replicate(1, V_2.rows())).transpose();
		return l;
	}

	double icp_step(Eigen::MatrixXd const & V_1, Eigen::MatrixXd & V_2, Eigen::MatrixXd const & N_1)
	{
		// Actual icp step, with normal

		// subsample
		Eigen::MatrixXd V_subset;
		Eigen::MatrixXd N_subset;
		if (prob < 1)
		{
			std::vector<int> indices;
			subsample(V_1.rows(), prob, indices);
			V_subset.resize(indices.size(), 3);
			N_subset.resize(indices.size(), 3);

			for (int i = 0; i < indices.size(); ++i)
			{
				V_subset.row(i) = V_1.row(indices[i]);
				N_subset.row(i) = N_1.row(indices[i]);
			}
		} else
		{
			V_subset = V_1;
			N_subset = N_1;
		}
		
		// setup coorespondance
		Eigen::MatrixXd V_ref;
		Eigen::MatrixXd V_NP;
		Eigen::MatrixXd N_ref;
		std::vector<int>  selection_ref;
		register_closest_point(V_subset, V_2, N_subset, V_ref, V_NP, N_ref);
		
		// calculate transformation
		Eigen::Matrix3d R;
		Eigen::Vector3d t;
		double l = estimate_rotation_translation(V_ref, N_ref, V_NP, R, t);
		if (mode != 2)
		{
			loss.push_back(l);
		}
		V_2 = (R*(V_2.transpose()) + t.replicate(1, V_2.rows())).transpose();
		return l;
	}

	void icp_step()
	{
		// public handle
		// check for mode and do one step accordingly
		if (mode != 2)
		{
			if (use_normal)
			{
				calculate_vertex_normal_flann(V_set_1, N_set_1);
				icp_step(V_set_1, V_set_2, N_set_1);
				return;
			}
			icp_step(V_set_1, V_set_2);
			return;
		}
		double l = 0;
		double li;
		if (use_normal)
		{
			for (int k = 1; k < 5; ++k)
			{
				std::cout << "Calculating normal for mesh " << k << "/4" << std::endl;
				calculate_vertex_normal_flann(V_set_vec[k], N_set_vec[k]);
			}
		}
		for (int i = 2; i > 0; --i)
		{
			if (use_normal)
			{
				li = icp_step(V_set_vec[i-1], V_set_vec[i], N_set_vec[i-1]);
			} else{
				li = icp_step(V_set_vec[i-1], V_set_vec[i]);
			}
			l += li;
			que5_loss[i-1].push_back(li);
		}
		if (use_normal)
		{
			li = icp_step(V_set_vec[4], V_set_vec[3], N_set_vec[4]);
		}else {
			li = icp_step(V_set_vec[4], V_set_vec[3]);
		}
		que5_loss[2].push_back(li);
		l += li;
		if (use_normal)
		{
			li = icp_step(V_set_vec[0], V_set_vec[4], N_set_vec[0]);
		} else
		{
			li = icp_step(V_set_vec[0], V_set_vec[4]);
		}
		que5_loss[3].push_back(li);
		l += li;
		std::cout << "Total loss = " << l << std::endl;
		loss.push_back(l);
	}

	void output_loss()
	{
		if (loss.size() == 0)
		{
			return;
		}
		std::ofstream outputFile("../output/loss_mode_"+std::to_string(mode)+".csv");
		for (int i = 0; i < loss.size(); ++i)
		{
			outputFile << loss[i] << ", ";
		}
		if (mode == 2)
		{
			for (int i = 0; i < 4; ++i)
			{
				outputFile << "\n";
				for (int j = 0; j < que5_loss[i].size(); ++j)
				{
					outputFile << que5_loss[i][j] << ", ";
				}
			}
		}
		outputFile.close();
	}

private:

};

MyContext g_myctx;


bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{

	std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
	if (key=='q' || key=='Q')
	{
		exit(0);
	}
	return false;
}


int main(int argc, char *argv[])
{
	// Init the viewer
	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	// menu variable Shared between two menus
	double doubleVariable = 0.1f; 

	// Add content to the default menu window via defining a Lambda expression with captures by reference([&])
	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::InputDouble("double", &doubleVariable, 0, 0, "%.4f");

			// ... or using a custom callback
			static bool boolVariable = true;
			if (ImGui::Checkbox("bool", &boolVariable))
			{
				// do something
				std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
			}

			// Expose an enumeration type
			enum Orientation { Up = 0, Down, Left, Right };
			static Orientation dir = Up;
			ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");

			// We can also use a std::vector<std::string> defined dynamically
			static int num_choices = 3;
			static std::vector<std::string> choices;
			static int idx_choice = 0;
			if (ImGui::InputInt("Num letters", &num_choices))
			{
				num_choices = std::max(1, std::min(26, num_choices));
			}
			if (num_choices != (int)choices.size())
			{
				choices.resize(num_choices);
				for (int i = 0; i < num_choices; ++i)
					choices[i] = std::string(1, 'A' + i);
				if (idx_choice >= num_choices)
					idx_choice = num_choices - 1;
			}
			ImGui::Combo("Letter", &idx_choice, choices);

			// Add a button
			if (ImGui::Button("Print Hello", ImVec2(-1, 0)))
			{
				std::cout << "Hello\n";
			}
		}
	};

	// Add additional windows via defining a Lambda expression with captures by reference([&])
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(300, 300), ImGuiSetCond_FirstUseEver);
		ImGui::Begin( "MyProperties", nullptr, ImGuiWindowFlags_NoSavedSettings );

		// point size
		// [event handle] if value changed
		if (ImGui::InputFloat("point_size", &g_myctx.point_size))
		{
			std::cout << "point_size changed\n";
			viewer.data().point_size = g_myctx.point_size;
		}

		// line width
		// [event handle] if value changed
		if(ImGui::InputFloat("line_width", &g_myctx.line_width))
		{
			std::cout << "line_width changed\n";
			viewer.data().line_width = g_myctx.line_width;
		}

		// mode
		// [event handle] if value changed
		if (ImGui::SliderInt("mode", &g_myctx.mode, 0,2))
		{
			if (g_myctx.mode == 0)
			{
				g_myctx.rotation = 45;
				g_myctx.noise_std =  0;
			} else if (g_myctx.mode == 1)
			{
				g_myctx.rotation = 5;
				g_myctx.noise_std =  0;
			} else if (g_myctx.mode == 2)
			{
				g_myctx.rotation = 0;
				g_myctx.noise_std =  0;
			}
			g_myctx.init_display(viewer);
		}

		//mode-text
		if (g_myctx.mode==0)
		{
			ImGui::Text("mode: ICP demo");
		} else if (g_myctx.mode == 1) {
			ImGui::Text("mode: ICP Same mesh rotated ");
		} else if (g_myctx.mode == 2){
			ImGui::Text("mode: ICP five meshs");
		}

		// subsampling rate
		ImGui::SliderFloat("Subsampling rate", &g_myctx.prob, 0.1, 1, "%.1f");

		// noise to source
		if (ImGui::InputInt("x0.0001 Noise std", &g_myctx.noise_std))
		{
			g_myctx.init_display(viewer);
		}

		// initial rotation
		if (ImGui::InputInt("Initial Rotation", &g_myctx.rotation))
		{
			g_myctx.init_display(viewer);
		}

		// initial translation
		if (ImGui::Checkbox("Initial Translation", &g_myctx.init_trans))
		{
			g_myctx.init_display(viewer);
		}

		// use normal
		ImGui::Checkbox("Use normal in ICP", &g_myctx.use_normal);

		// use normal
		if (ImGui::Checkbox("Shade normal", &g_myctx.shade))
		{
			g_myctx.reset_display(viewer);
		}

		// Do one iteration of ICP
		if (ImGui::Button("Step", ImVec2(-1, 0)))
		{
			g_myctx.icp_step();
			g_myctx.reset_display(viewer);
		}

		// Do ten iteration of ICP
		if (ImGui::Button("10 Steps", ImVec2(-1,0)))
		{
			for (int i = 0; i < 10; ++i) {
				std::cout << "Doing step " << i+1 <<"/10" << std::endl;
				g_myctx.icp_step();
			}
			g_myctx.reset_display(viewer);
		}

		// Output loss of each iteration to file
		if (ImGui::Button("Ouput Loss", ImVec2(-1,0)))
		{
			g_myctx.output_loss();
		}

		ImGui::End();
	};


	// registered a event handler
	viewer.callback_key_down = &key_down;

	g_myctx.init_display(viewer);

	// Call GUI
	viewer.launch();

}
