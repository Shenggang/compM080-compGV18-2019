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

	MyContext() :notexture(true),mode(0),render_index(0),rank(5),implicit(false),step_size(0.0000001),noise_level(0)
	{

	}
	~MyContext() {}

	Eigen::MatrixXd V_original;
	Eigen::MatrixXd m_V;
	Eigen::MatrixXi m_F;
	Eigen::MatrixXd m_C;
	Eigen::MatrixXd eigen_mesh;
	Eigen::SparseMatrix<double> Area;

	int const max_rank = 50;

	char const *directory[2] = {"../example_meshes/", "../example_meshes_notexture/"};
	char const *meshname;
	char const *prev_mesh;
	bool notexture;

	bool computed = false;
	int mode;
	int render_index;
	int rank;

	bool implicit;
	double step_size;
	double noise_level;
	

	void concate(Eigen::MatrixXd const & VA, Eigen::MatrixXi const & FA, Eigen::MatrixXd const & VB, Eigen::MatrixXi const & FB,
		Eigen::MatrixXd & out_V, Eigen::MatrixXi & out_F)
	{

		out_V.resize(VA.rows() + VB.rows(), VA.cols());
		out_V << VA, VB;
		out_F.resize(FA.rows() + FB.rows(), FA.cols());
		out_F << FA, (FB.array() + VA.rows());
		
	}

	void smooth()
	{

	}

	void compute_color()
	{
		Eigen::SparseMatrix<double> L, Area, AreaInv;
		Eigen::VectorXd result;
		compute_area(m_V, m_F, Area);

		AreaInv.resize(m_V.rows(), m_V.rows());
		AreaInv.setZero();
		for (int i = 0; i < Area.rows(); ++i)
		{
			AreaInv.insert(i,i) = 1/Area.coeff(i,i);
		}

		if (render_index == 0)
		{
			compute_gaussian_curvature(m_V, m_F, result);
			//result = AreaInv * result;
		} else
		{
			if (render_index == 1)
				compute_uniform_matrix(m_V,m_F,L);
			else
			{
				compute_cotan_matrix(m_V,m_F,L);
				L = AreaInv*L;
			}
			compute_H(m_V, L, result);
		}
		result = result.array() *100/(result.maxCoeff() - result.minCoeff());
		result = result.array() - result.minCoeff();
		//replace by color scheme
		igl::parula(result, false, m_C);
	}

	void reset_display(igl::opengl::glfw::Viewer& viewer)
	{
		viewer.data().clear(); 
		// hide default wireframe
		viewer.data().show_lines = 0;
		viewer.data().show_overlay_depth = 1; 

		// add mesh
		viewer.data().set_mesh(m_V, m_F);
		viewer.core.align_camera_center(m_V, m_F); 
	}

	void reset_display(igl::opengl::glfw::Viewer& viewer, Eigen::MatrixXd const &m_C)
	{
		reset_display(viewer);
		viewer.data().set_colors(m_C);
	}

	void reload_mesh(igl::opengl::glfw::Viewer& viewer)
	{
		std::string fullDirectory = std::string(directory[notexture])+std::string(meshname);
		igl::readOBJ(fullDirectory, V_original, m_F);
		m_V = V_original;
		reset_display(viewer);
	}

	void render(igl::opengl::glfw::Viewer &viewer)
	{
		compute_color();
		reset_display(viewer, m_C);
	}

	void reconstruct(igl::opengl::glfw::Viewer &viewer)
	{
		if (!computed)
		{
			// compute first 100 eigenvectors
			std::cout << "Computing eigen decompositions" <<std::endl;
			compute_eigendecomposition(V_original, m_F, max_rank, eigen_mesh, Area);
			computed = true;
		}

		// reconstruct mesh
		std::cout << "Reconstructing mesh" << std::endl;
		m_V.setZero();
		for (int i = 0; i < rank-1; ++i)
		{
			int index = max_rank - i - 1;
			for (int j = 0; j < 3; ++j)
			{
				m_V.col(j) += V_original.col(j).dot(Area*eigen_mesh.col(index)) * eigen_mesh.col(index);
			}
		}
		reset_display(viewer);
	}

	void smooth(igl::opengl::glfw::Viewer &viewer)
	{
		Eigen::SparseMatrix<double> M,C;
		compute_area(m_V, m_F, M);
		compute_cotan_matrix(m_V, m_F, C);

		if (implicit)
		{
			//implicit smoothing
			Eigen::SparseMatrix<double> A;
			A = M - step_size*C;
			Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
			//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
			solver.analyzePattern(A);
			solver.factorize(A);
			if (solver.info() != Eigen::Success) { std::cout << "decomposition failed\n";exit(1); }
			Eigen::MatrixXd b;
			b = M*m_V;
			for (int i = 0; i < 3; i++)
			{
				// solve Ax=b 
				// ...
				Eigen::VectorXd solution = solver.solve(b.col(i));
				if (solver.info() != Eigen::Success) {std::cout << "solving failed\n";exit(1);}
				m_V.col(i) = solution;
			}
		} else
		{
			//explicit smoothing
			M = M.cwiseInverse();
			m_V = m_V + step_size*M*C*m_V;
		}
		reset_display(viewer);
	}

	void add_noise(igl::opengl::glfw::Viewer &viewer)
	{
		std::default_random_engine generator;
  		std::normal_distribution<double> distribution(0, noise_level);
		Eigen::MatrixXd temp;
		temp.resize(V_original.rows(), V_original.cols());
		for (int i = 0; i < m_V.rows(); ++i)
		{
			for (int j = 0; j < m_V.cols(); ++j)
			{
				temp(i,j) = V_original(i,j) + distribution(generator);
			}
		}
		m_V = temp;
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
	char const *filenames[8] = {"bunny.obj", "camel.obj", "camelhead.obj", "cow_manifold.obj", "cow_small_manifold.obj", "cow.obj", "cube.obj", "dragon.obj"};
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
				g_myctx.computed = false;
				g_myctx.prev_mesh = g_myctx.meshname;
			}
		}

		if (ImGui::Checkbox("No Texture", &g_myctx.notexture))
		{
			g_myctx.reload_mesh(viewer);
		}

		// mode selection
		ImGui::SliderInt("mode", &g_myctx.mode, 0,2);

		if (g_myctx.mode == 0)
		{
			//mode 0, show curvatures
			const char* listbox_items[] = { "Gaussian Curvature" , "Mean Curvature (Uniform)", "Mean Curvature (Laplace-Beltrami)" };

			ImGui::ListBox("listbox\n(single select)", &g_myctx.render_index, listbox_items, IM_ARRAYSIZE(listbox_items), 4);

			if (ImGui::Button("Render"))
			{
				g_myctx.render(viewer);
			}
		} else if (g_myctx.mode == 1)
		{
			//mode 1, Low-rank approximation
			ImGui::Text("Low-rank approximation");
			ImGui::SliderInt("Reconstruction Rank", &g_myctx.rank, 1, g_myctx.max_rank);
			
			if (ImGui::Button("Reconstruct"))
			{
				g_myctx.reconstruct(viewer);
			}
		} else if (g_myctx.mode == 2)
		{
			//mode 2, mesh smoothing
			ImGui::Text("Mesh smoothingn");
			ImGui::Checkbox("Implicit smoothing", &g_myctx.implicit);
			ImGui::InputDouble("Step Size", &g_myctx.step_size, 0.0000001, 0.000001, "%.7f");
			if (ImGui::Button("Smooth"))
			{
				g_myctx.smooth(viewer);
			}
			ImGui::InputDouble("Noise level", &g_myctx.noise_level, 0.00001, 0.0001, "%.5f");
			if (ImGui::Button("Add noise"))
			{
				g_myctx.add_noise(viewer);
			}
		}
		ImGui::End();
	};


	// registered a event handler
	viewer.callback_key_down = &key_down;

	g_myctx.reload_mesh(viewer);

	// Call GUI
	viewer.launch();

}
