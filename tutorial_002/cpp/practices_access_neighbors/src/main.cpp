#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/vertex_triangle_adjacency.h>
#include <imgui/imgui.h>
#include <iostream>
#include "mytools.h"



class MyContext
{
public:

	MyContext() :nv_len(0), point_size(20), line_width(5), sel_vidx(7), mode(0)
	{
	
		//initial vertices and faces
		get_cube(m_V, m_F);

		m_num_vex = m_V.rows();

		// calculate face normal
		igl::per_face_normals(m_V, m_F, m_FN);

		// build adjacent matrix 
		// Inputs:
		//  number of vertices
		//  face def.
		//
		// Outputs:
		//   VF   #V list of lists of incident faces (adjacency list) 
		//   VFi  #V list of lists    of  index of incidence     within incident faces listed
		//        vector<vector<int>>        VI                    triangle_face 
		std::vector<std::vector<int> > VF;
		std::vector<std::vector<int> > VFi;
		igl::vertex_triangle_adjacency(m_V.rows(), m_F, VF, VFi);

		// initial fake vertex normal
		m_fake_VN.resize(m_V.rows(),m_V.cols());
		for (int i=0 ; i < m_V.rows(); i++)
		{
			m_fake_VN.row(i) = m_FN.row(VF[i][0]);
		}
		
		// calculate VN
		calculate_vertex_normal(m_V, m_F, m_FN, m_smooth_VN);


		// save 
		m_VF = std::move(VF);
		m_VFi = std::move(VFi);

	}
	~MyContext() {}

	Eigen::MatrixXd m_V;
	Eigen::MatrixXi m_F;
	Eigen::Matrix<double, Eigen::Dynamic, 3> m_FN;
	
	Eigen::MatrixXd m_fake_VN;
	Eigen::MatrixXd m_smooth_VN;

	std::vector<std::vector<int> > m_VF;
	std::vector<std::vector<int> > m_VFi;

	int m_num_vex;
	float nv_len;
	float point_size;
	float line_width;
	int sel_vidx;
	int mode;

	void reset_display(igl::opengl::glfw::Viewer& viewer)
	{
		
		viewer.data().clear();
		viewer.data().set_mesh(m_V, m_F);
		viewer.core.align_camera_center(m_V, m_F);

		// hide default wireframe
		viewer.data().show_lines = 1;
		viewer.data().show_overlay_depth = 0;
		viewer.data().show_faces = 0;

		//======================================================================
		// visualize adjacent vertices
		{
			Eigen::MatrixXd adj_vex;
			adj_vex.resize(m_VF[sel_vidx].size() * 2, 3);

			int count = 0;
			// get adjacent faces & vertices
			for (size_t i = 0; i < m_VF[sel_vidx].size(); i++)
			{
				int face_idx = m_VF[sel_vidx][i];
				int v_local_idx = m_VFi[sel_vidx][i];

				for (int iv = 0; iv < 3; iv++)
				{
					if (iv != v_local_idx)
					{
						Eigen::RowVector3d const vex = m_V.row(m_F(face_idx, iv));
						adj_vex.row(count) = vex;
						count++;
					}
				}
			} 
			
			//mark points
			viewer.data().add_points(m_V.row(sel_vidx), Eigen::RowVector3d(1, 0, 0));
			viewer.data().add_points(adj_vex, Eigen::RowVector3d(1, 1, 0));

			// add links
			Eigen::MatrixXd EV1 = adj_vex;
			Eigen::MatrixXd EV2(adj_vex.rows(), 3);
			EV2.setZero();
			EV2 = EV2.rowwise() + m_V.row(sel_vidx);

			viewer.data().add_edges(EV1, EV2, Eigen::RowVector3d(0, 0, 1));
		} 
		
		//======================================================================
		// add normal lines
		{
			Eigen::MatrixXd EV1(m_V);
			Eigen::MatrixXd EV2;
			if (mode==0)
			{
				EV2 = m_V + m_fake_VN * nv_len;
			}
			else {
				// show real VN
				EV2 = m_V + m_smooth_VN * nv_len;

			}

			viewer.data().add_edges(EV1, EV2, Eigen::RowVector3d(1, 1, 1));

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
		ImGui::SetNextWindowSize(ImVec2(300, 160), ImGuiSetCond_FirstUseEver);
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
		if (ImGui::SliderInt("mode", &g_myctx.mode, 0, 1))
		{
			g_myctx.reset_display(viewer);
		}

		//mode-text
		if (g_myctx.mode==0)
		{
			ImGui::Text("fake vertex normal");
		}
		else {
			ImGui::Text("vertex normal");
		}

		ImGui::End();
	};


	// registered a event handler
	viewer.callback_key_down = &key_down;

	g_myctx.reset_display(viewer);

	// Call GUI
	viewer.launch();

}
