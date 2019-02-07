#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/vertex_triangle_adjacency.h>
#include <imgui/imgui.h>
#include <igl/readPLY.h>
#include <igl/readOff.h>
#include <igl/file_exists.h>
#include <iostream>
#include "mytools.h"



class MyContext
{
public:

	MyContext() :nv_len(5), point_size(5), line_width(1), mode(0), slice_y(0.5), energy_color(5)
	{

	}
	~MyContext() {}

	Eigen::MatrixXd m_V;
	Eigen::MatrixXi m_F;	 
	Eigen::MatrixXd m_VN;
	Eigen::MatrixXd input_pts;
	Eigen::MatrixXd input_color;
	Eigen::VectorXd sdf_param;

	int m_num_vex;
	float nv_len;
	float point_size;
	float line_width;

	Eigen::Vector3d max_pts;
	Eigen::Vector3d min_pts; 
	float slice_y; 
	float energy_color;
 	int mode;


	void build_slice_and_samples(Eigen::MatrixXd & out_splane_V, Eigen::MatrixXi & out_splane_F,
		Eigen::MatrixXd & out_samples )
	{
		Eigen::MatrixXd splane_v(4, 3);
		/* 
		       x
			0------3
			| \    |
			|   \  | z
			|     \|
		    1------2
		*/

		splane_v.rowwise() = min_pts.transpose();
		splane_v.row(2) = max_pts;

		splane_v(1,2) = max_pts(2);
		splane_v(3,0) = max_pts(0);

		//sey Y
		float y_val = (1 - slice_y)*max_pts.y() + slice_y * min_pts.y();
		splane_v.col(1).setConstant(y_val);
		 
		//---------------------------------------------------------
		// build plane's mesh
		Eigen::MatrixXi splane_F(2, 3);
		splane_F << 0, 1, 2,
					0, 2, 3;

		out_splane_V = splane_v;
		out_splane_F = splane_F;

		//---------------------------------------------------------

		Eigen::Vector3d pl_v0 = splane_v.row(0);
		Eigen::Vector3d pl_v2 = splane_v.row(2);

		//250*250
		float step_size = (pl_v2.x() - pl_v0.x()) / 100;

		size_t num_x = (pl_v2.x() - pl_v0.x()) / step_size +1;
		size_t num_z = (pl_v2.z() - pl_v0.z()) / step_size +1;

		out_samples.resize(num_x*num_z, 3);


		for (int xi = 0; xi < num_x; xi++)
		{
			for (int zi = 0; zi < num_z; zi++)
			{
				// xi,zi
				out_samples.row(xi*num_z +zi)<< pl_v0.x()+xi* step_size,y_val, pl_v0.z()+zi* step_size;

			}
		}


	}



	void concate(Eigen::MatrixXd const & VA, Eigen::MatrixXi const & FA, Eigen::MatrixXd const & VB, Eigen::MatrixXi const & FB,
		Eigen::MatrixXd & out_V, Eigen::MatrixXi & out_F	)
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

		//======================================================================
		// build  a slice plane + sample points
		Eigen::MatrixXd splane_V;
		Eigen::MatrixXi splane_F;

		Eigen::MatrixXd splane_spts;
		build_slice_and_samples(splane_V, splane_F, splane_spts);


		//======================================================================

		viewer.data().line_width = line_width;
		viewer.data().point_size = point_size;

		if (mode == 0 )
		{
			viewer.data().set_points(input_pts, input_color);

			viewer.data().set_mesh(m_V, m_F);
			viewer.core.align_camera_center(m_V, m_F);

			// add normal lines
			{
				Eigen::MatrixXd EV1(m_V);
				Eigen::MatrixXd EV2;
				if (mode == 0)
				{
					EV2 = m_V + m_VN * nv_len;
				}
				else {
					// show real VN
					EV2 = m_V + m_VN * nv_len;
				}

				viewer.data().add_edges(EV1, EV2, Eigen::RowVector3d(1, 1, 1));

			}

		}
		else if (mode == 1)
		{
			// show two meshes

			// Concatenate (VA,FA) and (VB,FB) into (V,F)
			
			Eigen::MatrixXd dis_V;
			Eigen::MatrixXi dis_F;
			concate(m_V, m_F, splane_V, splane_F, dis_V, dis_F);

			viewer.data().set_mesh(dis_V, dis_F);
			viewer.core.align_camera_center(dis_V, dis_F);

		}
		else if (mode == 2) { 
			// calculate SDF
			Eigen::VectorXd sdf_vals;
			calculate_SDF(input_pts, sdf_param, splane_spts, sdf_vals);

			//change sdf_sspts to vertex color 
			Eigen::MatrixXd sdf_color(sdf_vals.rows(), 3);
			sdf_color.setZero();

			const Eigen::VectorXd ct_1 = Eigen::VectorXd::Constant(sdf_color.rows(), 1);

			// set to white(1,1,1) if input_sdist[i] == 0 else remain the same
			sdf_color.col(0) = (sdf_vals.array() == 0).select(ct_1, sdf_color.col(0));
			sdf_color.col(1) = (sdf_vals.array() == 0).select(ct_1, sdf_color.col(1));
			sdf_color.col(2) = (sdf_vals.array() == 0).select(ct_1, sdf_color.col(2));

			Eigen::VectorXd  energy = sdf_vals.array().abs();
			energy = energy.array() / (energy_color);
			energy = (energy.array() > 1).select(ct_1, energy);
			energy = 1 - energy.array();

			//std::cout << "max energy=" << sdf_vals.array().maxCoeff();
			//std::cout << "max energy=" << sdf_vals.array().minCoeff();

			// set to red(e,1,0) if input_sdist[i] < 0 else remain the same
			sdf_color.col(0) = (sdf_vals.array() < 0).select(energy, sdf_color.col(0));

			// set to green(0,e,0) if input_sdist[i] > 0 else remain the same
			sdf_color.col(1) = (sdf_vals.array() > 0).select(energy, sdf_color.col(1));


			viewer.data().add_points(splane_spts, sdf_color);

			viewer.data().set_mesh(m_V, m_F);
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

void get_example_mesh(Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::MatrixXd & VN)
{

	std::vector<const char *> cands{ "../../../data/example_small.off", "../../../../data/example_small.off", "../../data/example_small.off" };

	bool found = false;
	for (const auto & val : cands)
	{
		if ( igl::file_exists(val) )
		{	
			std::cout << "loading example mesh from:" << val << "\n";

			if (igl::readOFF(val, V,F)) {
				igl::per_vertex_normals(V, F, VN);
				found = 1;
				break;
			}
		}
	}

	if (!found) {
		std::cout << "cannot locate ../../../data/example_small.off\n";
		exit(1);
	}
	

}





int main(int argc, char *argv[])
{
	//------------------------------------------
	// load data at...
	//   tutorial_003\data
	Eigen::MatrixXd V;
	Eigen::MatrixXd VN;
	Eigen::MatrixXi F;  
	get_example_mesh(V, F, VN); 

	Eigen::MatrixXd input_pts;
	Eigen::VectorXd input_sdist;
	create_input_pcd(V,VN, input_pts,input_sdist);

	// extract SDF
	Eigen::VectorXd lambda_c;
	fit_SDF(input_pts,input_sdist, lambda_c);


	//------------------------------------------
	// for visualization

	g_myctx.m_V = V;
	g_myctx.m_F = F;
	g_myctx.m_VN = VN;
	
	g_myctx.max_pts = input_pts.colwise().maxCoeff();
	g_myctx.min_pts = input_pts.colwise().minCoeff();
	
	g_myctx.slice_y = 0.5;  

	g_myctx.input_pts = input_pts;
	g_myctx.input_color.resizeLike(input_pts);
	g_myctx.input_color.setZero();

	g_myctx.sdf_param = lambda_c;

	const Eigen::VectorXd ct_1 = Eigen::VectorXd::Constant(input_sdist.rows(), 1);
	const Eigen::VectorXd ct_0 = Eigen::VectorXd::Constant(input_sdist.rows(), 0);
	
	// set to white(1,1,1) if input_sdist[i] == 0 else remain the same
	g_myctx.input_color.col(0) = (input_sdist.array() == 0).select(ct_1, g_myctx.input_color.col(0));
	g_myctx.input_color.col(1) = (input_sdist.array() == 0).select(ct_1, g_myctx.input_color.col(1));
	g_myctx.input_color.col(2) = (input_sdist.array() == 0).select(ct_1, g_myctx.input_color.col(2));

	// set to red(1,0,0) if input_sdist[i] < 0 else remain the same
	g_myctx.input_color.col(0) = (input_sdist.array() < 0).select(ct_1, g_myctx.input_color.col(0));

	// set to red(0,1,0) if input_sdist[i] > 0 else remain the same
	g_myctx.input_color.col(1) = (input_sdist.array() > 0).select(ct_1, g_myctx.input_color.col(1));


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

		// length of normal line
		// [event handle] if value changed
		//if (ImGui::InputFloat("nv_length", &g_myctx.nv_len))
		if (ImGui::SliderFloat("nv_length", &g_myctx.nv_len, 0, 50, "%.1f"))
		{
			//pass
			g_myctx.reset_display(viewer);
		}


		//mode - List box
		const char* listbox_items[] = { "Mesh", "Slice Plane", "SDF"};

		if (ImGui::ListBox("listbox\n(single select)", &g_myctx.mode, listbox_items, IM_ARRAYSIZE(listbox_items), 4))
		{

			g_myctx.reset_display(viewer);
		}

		if (ImGui::SliderFloat("Y_slice_%", &g_myctx.slice_y, 0, 1, "%.2f"))
		{
			//pass
			g_myctx.reset_display(viewer);
		}
		if (ImGui::SliderFloat("energy_color", &g_myctx.energy_color, 0.01, 10, "%.2f"))
		{
			//pass
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
