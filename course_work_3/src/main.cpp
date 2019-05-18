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
#include "seamster.h"



class MyContext
{
public:

	MyContext() :notexture(true),line_width(1),point_size(5),mode(0)
	{

	}
	~MyContext() {}

	Eigen::MatrixXd V_original;
	Eigen::MatrixXd m_V, V_cut, V_flat_uniform, V_flat_LB, V_flat_mean, EV1, EV2;
	Eigen::MatrixXi m_F, F_cut;
	Eigen::MatrixXd m_C;
	double distortion[3];
	Eigen::SparseMatrix<double> Area;

	char const *directory[2] = {"../example_meshes/","../example_meshes_notexture/"};
	char const *meshname;
	char const *prev_mesh;
	bool notexture;
	bool computed = false;
	bool seam_cut = false;
	float line_width;
	float point_size;
	int mode;

	void concate(Eigen::MatrixXd const & VA, Eigen::MatrixXi const & FA, Eigen::MatrixXd const & VB, Eigen::MatrixXi const & FB,
		Eigen::MatrixXd & out_V, Eigen::MatrixXi & out_F)
	{
		out_V.resize(VA.rows() + VB.rows(), VA.cols());
		out_V << VA, VB;
		out_F.resize(FA.rows() + FB.rows(), FA.cols());
		out_F << FA, (FB.array() + VA.rows());
	}

	void reset_display(igl::opengl::glfw::Viewer& viewer)
	{
		viewer.data().clear(); 
		// hide default wireframe
		viewer.data().show_lines = 1;
		viewer.data().show_overlay_depth = 1; 

		viewer.data().line_width = line_width;
		viewer.data().point_size = point_size;
		
		// add mesh
		if (mode == 0)
		{
			viewer.data().set_mesh(m_V, m_F);
			viewer.core.align_camera_center(m_V, m_F); 
			if (seam_cut)
			{
				// show seam edges
				viewer.data().clear();
				viewer.data().set_mesh(V_cut, F_cut);
				viewer.core.align_camera_center(V_cut, F_cut); 
				viewer.data().add_edges(EV1, EV2, Eigen::RowVector3d(0, 0, 1));
			}
		}
		else if (mode == 1) 
		{
			viewer.data().set_mesh(V_flat_uniform, F_cut);
			viewer.core.align_camera_center(V_flat_uniform, F_cut);
		}
		else if (mode == 2)
		{
			viewer.data().set_mesh(V_flat_LB, F_cut);
			viewer.core.align_camera_center(V_flat_LB, F_cut); 
		}
		else if (mode == 3)
		{
			viewer.data().set_mesh(V_flat_mean, F_cut);
			viewer.core.align_camera_center(V_flat_mean, F_cut); 
		}
	}

	void reset_display(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd const &m_C)
	{
		reset_display(viewer);
		viewer.data().set_colors(m_C);
	}

	void reload_mesh(igl::opengl::glfw::Viewer& viewer)
	{
		std::string fullDirectory = std::string(directory[0])+std::string(meshname);
		igl::readOBJ(fullDirectory, V_original, m_F);
		m_V = V_original;
		mode = 0;
		computed = false;
		seam_cut = false;
		V_cut.setZero();
		F_cut.setZero();
		V_flat_LB.setZero();
		V_flat_uniform.setZero();
		reset_display(viewer);
	}

	void reset(igl::opengl::glfw::Viewer &viewer)
	{
		m_V = V_original;
		reset_display(viewer);
	}

	void flatten(igl::opengl::glfw::Viewer& viewer)
	{
		Eigen::MatrixXd FA_cut, FA_uni, FA_LB, FA_mean;
		compute_face_angle(V_cut, F_cut, FA_cut);
		// uniform laplace
		tutte_flattening_uniform(V_cut, F_cut, V_flat_uniform);
		compute_face_angle(V_flat_uniform, F_cut, FA_uni);
		// Laplace-Beltrami
		tutte_flattening_LB(V_cut, F_cut, V_flat_LB);
		compute_face_angle(V_flat_LB, F_cut, FA_LB);
		// Mean laplce
		tutte_flattening_mean(V_cut, F_cut, V_flat_mean);
		compute_face_angle(V_flat_mean, F_cut, FA_mean);

		// Compute distortion
		distortion[0] = compute_distortion(FA_cut, FA_uni);
		distortion[1] = compute_distortion(FA_cut, FA_LB);
		distortion[2] = compute_distortion(FA_cut, FA_mean);

		computed = true;
		mode = 1;
		reset_display(viewer);
	}

	void cut(igl::opengl::glfw::Viewer& viewer)
	{

		// cut mesh
		EdgeTopology et(m_V, m_F);
		// compute seam edges
		vector<int>  Terminals(0);
		find_seam_terminals(m_V, m_F, Terminals);
		cout << Terminals.size() << endl;
		Eigen::VectorXi b_bex_index;
		igl::boundary_loop(m_F, b_bex_index);
		vector<int> SEi;
		SEi.clear();
		if (b_bex_index.rows() <10)
		{
			cout << "Compute seam" << endl;
			compute_seam(m_V, m_F, et, Terminals, SEi);
		}
		cout << "Seam found" << endl;
		// set seam edges
		EV1.resize(SEi.size(),3);
		EV2.resize(SEi.size(),3);
		for (int i = 0; i < SEi.size(); ++i)
		{
			EV1.row(i) = m_V.row(et.EV(SEi[i],0));
			EV2.row(i) = m_V.row(et.EV(SEi[i],1));
		}
		cut_mesh(m_V, m_F, et, SEi, V_cut, F_cut);

		seam_cut = true;
		mode = 0;
		reset_display(viewer);
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
	char const *filenames[4] = {"bunny.obj", "camel.obj", "camelhead.obj", "CrumpledDevelopable.obj"};
	g_myctx.meshname = filenames[0];
	g_myctx.prev_mesh = filenames[0];

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

		// mesh selection
		if (ImGui::BeginCombo("##combo", g_myctx.meshname))
		{
			for (int n = 0; n < IM_ARRAYSIZE(filenames); ++n)
			{
				bool is_selected = (g_myctx.meshname == filenames[n]);
				if (ImGui::Selectable(filenames[n], is_selected))
					g_myctx.meshname = filenames[n];
				if (is_selected)
					ImGui::SetItemDefaultFocus();
			}
			ImGui::EndCombo();
			if (g_myctx.prev_mesh != g_myctx.meshname)
			{
				g_myctx.reload_mesh(viewer);
				g_myctx.prev_mesh = g_myctx.meshname;
			}
		}

		if (g_myctx.computed)
		{
			if (ImGui::SliderInt("Mode", &g_myctx.mode, 0,3))
			{
				g_myctx.reset_display(viewer);
			}
			if (g_myctx.mode > 0)
			{
				ImGui::Text("Energy = %.4f", g_myctx.distortion[g_myctx.mode-1]);
			}
		}
		if (!g_myctx.seam_cut)
		{
			if (ImGui::Button("Cut"))
			{
				g_myctx.cut(viewer);
			}
		}

		if (!g_myctx.computed && g_myctx.seam_cut)
		{
			if (ImGui::Button("Flatten"))
			{
				g_myctx.flatten(viewer);
			}
		}

		// More functions
		ImGui::End();
	};


	// registered a event handler
	viewer.callback_key_down = &key_down;

	g_myctx.reload_mesh(viewer);

	// Call GUI
	viewer.launch();

}
