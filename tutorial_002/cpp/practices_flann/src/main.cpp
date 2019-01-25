#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/writeOFF.h>
#include <igl/file_exists.h>
#include <imgui/imgui.h>
#include <iostream>
#include <random>
#include "mytools.h"



class MyContext
{
public:

	MyContext() :nv_len(0.5), point_size(8), line_width(0.5), sel_vidx(7), mode(0)
	{
	
		//initial vertices and faces
		if (!igl::file_exists("../data/fandisk.off"))
		{
			std::cout << "[error] cannot locate model file at '../data/fandisk.off' \nPress any key to exit\n";
			char c;
			std::cin >> c;
			exit(1);
		}

		igl::readOFF("../data/fandisk.off", m_V, m_F);

		m_num_vex = m_V.rows();


		
		// calculate VN 
		igl::per_face_normals(m_V, m_F, m_FN);

		calculate_vertex_normal(m_V, m_F, m_FN, m_VN);

		calculate_vertex_normal_flann(m_V, m_F, m_VN_flann);


	}
	~MyContext() {}

	Eigen::MatrixXd m_V;
	Eigen::MatrixXi m_F; 

	Eigen::MatrixXd m_FN;
	Eigen::MatrixXd m_VN;
	Eigen::MatrixXd m_VN_flann;


	int m_num_vex;
	float nv_len;
	float point_size;
	float line_width;
	int sel_vidx;
	int mode;

	void reset_display(igl::opengl::glfw::Viewer& viewer)
	{

		static std::default_random_engine generator;
		
		viewer.data().clear();

		if (mode == 1 || mode == 2)
		{

			viewer.data().set_mesh(m_V, m_F);
			viewer.core.align_camera_center(m_V, m_F);
			// hide default wireframe
			viewer.data().show_lines = 0;

			//======================================================================
			// show normal lines
			Eigen::MatrixXd EV1(m_V);
			Eigen::MatrixXd EV2;

			if (mode == 1)
			{				
				EV2 = m_V + m_VN * nv_len;
			}
			else {

				EV2 = m_V + m_VN_flann * nv_len;
			}
			viewer.data().add_edges(EV1, EV2, Eigen::RowVector3d(1, 1, 1));

			viewer.core.align_camera_center(m_V,m_F);
		}
		else
		{
			// FLANN demo
			
			typedef nanoflann::KDTreeEigenMatrixAdaptor< Eigen::MatrixXd >  my_kd_tree_t;

			// build kd-tree uses vertices m_V
			// leaf_max_size:
			// "Large values mean that the tree will be built faster(since the tree will be smaller), but each query will be slower(since the linear search in the leaf is to be done over more points)."
			// "Small values will build the tree much slower(there will be many tree nodes), but queries will be faster... up to some point, since the "tree-part" of the search(logarithmic complexity) still has a significant cost."
			my_kd_tree_t mat_index( m_V, 50 /* max leaf */);
			mat_index.index->buildIndex();

			// select a random vertex
			std::uniform_int_distribution<int> distribution(0, m_V.rows());
			int rd_vidx = distribution(generator);

			Eigen::RowVector3d rd_pts = m_V.row(rd_vidx);

			// set K nearest samples
			const size_t par_K = 20;

			// create a query object
			std::vector<size_t> indexes(par_K);
			std::vector<double> dists_sqr(par_K);

			nanoflann::KNNResultSet<double> res(par_K);
			res.init(indexes.data(), dists_sqr.data());

			// find KNN
			//   SearchParams Note: The first argument (checks_IGNORED_) is ignored, but kept for compatibility
			//   checks_IGNORED_ was used to specify maximum # of leaf to be visited 
			mat_index.index->findNeighbors(res, rd_pts.data(), nanoflann::SearchParams(50));

			Eigen::MatrixXd nn_vex(indexes.size(),3);
			for (size_t i=0 ; i < indexes.size() ;i++)
			{
				nn_vex.row(i) = m_V.row(indexes[i]);
			}

			// show 
			viewer.data().add_points(rd_pts, Eigen::RowVector3d(0, 0, 1));
			viewer.data().add_points(nn_vex, Eigen::RowVector3d(1, 0, 0));

			Eigen::MatrixXd EV1(nn_vex);
			EV1.setZero();
			EV1 = EV1.rowwise() + rd_pts;
			viewer.data().add_edges(EV1, nn_vex, Eigen::RowVector3d(1, 1, 1));

			// add other points
			viewer.data().add_points(m_V, Eigen::RowVector3d(0, 0, 0));
			viewer.core.align_camera_center(nn_vex);
			viewer.data().show_overlay_depth = 1;
			viewer.data().show_overlay = 1;
		}


		//======================================================================

		viewer.data().line_width = line_width;
		viewer.data().point_size = point_size;

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
		ImGui::SetNextWindowSize(ImVec2(250, 300), ImGuiSetCond_FirstUseEver);
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

		// length of normal line
		// [event handle] if value changed
		//if (ImGui::InputFloat("nv_length", &g_myctx.nv_len))
		if (ImGui::SliderFloat("nv_length", &g_myctx.nv_len, 0, 1, "%.1f"))
		{
			//pass
			g_myctx.reset_display(viewer);
		}

		// vertex index
		if (ImGui::SliderInt("sel_vex_index", &g_myctx.sel_vidx, 0, g_myctx.m_num_vex-1))
		{
			g_myctx.reset_display(viewer);
		} 

		//mode
		if (ImGui::SliderInt("mode", &g_myctx.mode, 0,2))
		{
			g_myctx.reset_display(viewer);
		}

		//mode-text
		if (g_myctx.mode==0)
		{
			ImGui::Text("mode: KNN demo");
		}
		else if (g_myctx.mode == 1) { 
			ImGui::Text("mode: vertex normal ");
		}
		else {
			ImGui::Text("mode: vertex normal(flann)");
		}

		if (ImGui::Button("Random Query", ImVec2(-1, 0)))
		{
			g_myctx.reset_display(viewer);
		}
		ImGui::End();
	};


	// registered a event handler
	viewer.callback_key_down = &key_down;

	g_myctx.reset_display(viewer);

	// Call GUI
	viewer.launch();

}
