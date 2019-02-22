#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/vertex_triangle_adjacency.h>
#include <imgui/imgui.h>
#include <igl/readPLY.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/file_exists.h>
#include <igl/setdiff.h> 
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/LU>


#include <iostream>
#include "mytools.h"



class MyContext
{
public:

	MyContext() :nv_len(0), point_size(5), line_width(10), mode(0)
	{

	}
	~MyContext() {}

	Eigen::MatrixXd m_V;
	Eigen::MatrixXi m_F;	 
	Eigen::MatrixXd m_VN;
	Eigen::MatrixXd m_C;
	Eigen::MatrixXd input_pts;
	Eigen::MatrixXd bvex;
	Eigen::MatrixXi bedges; 
	Eigen::SparseMatrix<double> lmat;

	int m_num_vex;
	float nv_len;
	float point_size;
	float line_width;

 	int mode;

	void concate(Eigen::MatrixXd const & VA, Eigen::MatrixXi const & FA, Eigen::MatrixXd const & VB, Eigen::MatrixXi const & FB,
		Eigen::MatrixXd & out_V, Eigen::MatrixXi & out_F	)
	{

		out_V.resize(VA.rows() + VB.rows(), VA.cols());
		out_V << VA, VB;
		out_F.resize(FA.rows() + FB.rows(), FA.cols());
		out_F << FA, (FB.array() + VA.rows());
		
	}

	void smooth()
	{

	}

	void reset_display(igl::opengl::glfw::Viewer& viewer)
	{
		
		viewer.data().clear(); 
		// hide default wireframe
		viewer.data().show_lines = 0;
		viewer.data().show_overlay_depth = 1; 


		//======================================================================

		viewer.data().line_width = line_width;
		viewer.data().point_size = point_size;

		if (mode == 0 )
		{

			// add boundary vertices
			if (bvex.size()>0)
			{
				Eigen::MatrixXd vcolor;
				vcolor.resizeLike(bvex);
				vcolor.setZero();
				vcolor.col(0).setOnes();

				viewer.data().set_points(bvex, vcolor);

				// add edges
				const size_t NUM = bedges.rows();
				Eigen::MatrixXd EV1(NUM, 3);
				Eigen::MatrixXd EV2(NUM, 3);

				for (size_t i = 0; i < NUM; i++)
				{
					EV1.row(i) = m_V.row(bedges(i, 0));
					EV2.row(i) = m_V.row(bedges(i, 1));
				}

				viewer.data().add_edges(EV1, EV2, Eigen::RowVector3d(0, 0, 1));
			}



			// add mesh
			viewer.data().set_mesh(m_V, m_F);
			viewer.core.align_camera_center(m_V, m_F);


		}
		else if (mode == 1)
		{

			// add mesh
			viewer.data().set_mesh(m_V, m_F);
			viewer.data().set_colors(m_C);
			viewer.core.align_camera_center(m_V, m_F); 
		}


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

void get_example_mesh(std::string const meshname , Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::MatrixXd & VN)
{


	std::vector<const char *> cands{ 
		"../../data/", 
		"../../../data/",
		"../../../../data/",
		"../../../../../data/" };

	bool found = false;
	for (const auto & val : cands)
	{
		if ( igl::file_exists(val+ meshname) )
		{	
			std::cout << "loading example mesh from:" << val+ meshname << "\n";

			if (igl::readOFF(val+ meshname, V,F)) {
				igl::per_vertex_normals(V, F, VN);
				found = 1;
				break;
			}
			else {
				std::cout << "file loading failed " << cands[0] + meshname << "\n"; 
			}
		}
	}

	if (!found) {
		std::cout << "cannot locate "<<cands[0]+ meshname <<"\n";
		exit(1);
	}

}




int main(int argc, char *argv[])
{
	//------------------------------------------
	// load data  
	Eigen::MatrixXd V;
	Eigen::MatrixXd VN;
	Eigen::MatrixXi F;  

	get_example_mesh("camelhead.off", V, F, VN);
	//get_example_mesh("fertility.off", V, F, VN);
	//get_example_mesh("face_cut_1_asc.off", V, F, VN);

	//------------------------------------------
	// call your func.

	Eigen::MatrixXd bvex;
	get_boundary_vex(V, F, bvex);

	Eigen::MatrixXi bedges;
	get_boundary_edges(F, bedges);

	Eigen::VectorXd H;
	compute_H(V, F, H);

	//------------------------------------------
	// for visualization
	g_myctx.m_V = V;
	g_myctx.m_F = F;
	g_myctx.m_VN = VN;
	g_myctx.bedges = bedges;
	g_myctx.bvex = bvex;

	H = 100 * H.array() / (H.maxCoeff() - H.minCoeff());
	//replace by color scheme
	igl::parula(H, false, g_myctx.m_C);

	//------------------------------------------
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

		}
	};

	// Add additional windows via defining a Lambda expression with captures by reference([&])
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(250, 400), ImGuiSetCond_FirstUseEver);
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


		//mode - List box
		const char* listbox_items[] = { "Vex/Edge" , "H" };

		if (ImGui::ListBox("listbox\n(single select)", &g_myctx.mode, listbox_items, IM_ARRAYSIZE(listbox_items), 4))
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
