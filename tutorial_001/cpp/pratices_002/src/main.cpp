

#include "mytools.h"

#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
#include <igl/writePLY.h>
#include <igl/writeObj.h>
#include <igl/file_exists.h>

#include "sensorData.h"


class MyContext
{
public:

	MyContext():nv_len(1.0),mode(1),frame_idx(-1){}
	~MyContext() {}

	//display variables
	int mode;
	float nv_len;
	int frame_idx;
	float point_size;
	float line_width;

	//scan data
	std::vector<ml::mat4f> cam_poses; 
	std::vector<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> depths;
	std::vector<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> colors;
	std::vector<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> point_clouds;
	std::vector<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> point_clouds_rgb;
	Eigen::MatrixXf averaged_depth;
	Eigen::MatrixXf averaged_pcd;
	Eigen::MatrixXf averaged_pcd_rgb;

	void change_mode(igl::opengl::glfw::Viewer& viewer)
	{
		// call by event handler

		reset_display(viewer);

		std::cout << "change mode=" << mode <<"\n";
	}

	void set_camera(igl::opengl::glfw::Viewer & viewer)
	{
		Eigen::MatrixXd V = point_clouds[frame_idx].cast<double>();
		viewer.core.align_camera_center(V); 
	}

	void reset_display(igl::opengl::glfw::Viewer& viewer)
	{
		// reset display setting 

		viewer.data().point_size = point_size;
		viewer.data().line_width = line_width;

		// fix gui's camera 
		Eigen::Matrix3d mv;
		mv.setIdentity();
		mv.row(1) << 0, -1, 0;
		mv.row(2) << 0, 0, -1;

		if (mode == 1)
		{
			if (frame_idx!=-1 && frame_idx < point_clouds.size() ) 
			{
				Eigen::MatrixXd V = point_clouds[frame_idx].cast<double>();

				V = (mv * (V.transpose())).transpose();

				viewer.data().clear();
				viewer.data().add_points( V, point_clouds_rgb[frame_idx].cast<double>());

			}
		}
		else if (mode == 2)
		{
			Eigen::MatrixXd V = averaged_pcd.cast<double>();

			V = (mv * (V.transpose())).transpose(); 

			viewer.data().clear();
			viewer.data().add_points(V, averaged_pcd_rgb.cast<double>());

		}
	}

private:


};


MyContext g_myctx;

//=================================================================================


bool key_press(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
	std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;

	if (key == 'N' || key == 'n') {
		g_myctx.frame_idx++;

		if (g_myctx.frame_idx > g_myctx.point_clouds.size()) {
			g_myctx.frame_idx = 0;
		}

		g_myctx.reset_display(viewer);
	}

	return false;
}

void ShowHelpMarker(const char* desc)
{
	ImGui::TextDisabled("(?)");
	if (ImGui::IsItemHovered())
	{
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(desc);
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}
}

void ml2eigen(ml::mat4f const & mat, Eigen::Matrix4f & out_mat)
{
	//out_mat.resize();

	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			//row major
			out_mat(i, j) = mat.matrix2[i][j];
		}
	}

}

void unproj(unsigned int const WIDTH, unsigned int const HEIGHT,
	Eigen::MatrixXf const & depth_intr_inv,
	Eigen::MatrixXf const & color_extr_inv,
	Eigen::MatrixXf const & i_depth, Eigen::MatrixXf const & i_color,
	Eigen::MatrixXf & out_pts, Eigen::MatrixXf & out_pts_rgb)
{

	std::vector<Eigen::Vector3f> pts_list;
	std::vector<Eigen::Vector3f> pts_rgb_list;

	for (unsigned int i = 0; i < WIDTH*HEIGHT; i++) {

		float d = i_depth(i);

		Eigen::Vector3f rgb = i_color.row(i);

		//skip invalid depth 
		if (d != 0) {

			// row major 
			unsigned int x = i % WIDTH;
			unsigned int y = i / HEIGHT;

			Eigen::Vector4f dep_cam_pos = (depth_intr_inv * Eigen::Vector4f((float)x*d, (float)y*d, d, 0.0f));
			dep_cam_pos(3) = 1.0;

			// tx from depth camera space to color camera space 
			// it's the same in primesense sensor
			Eigen::Vector4f color_cam_pos = dep_cam_pos;
			color_cam_pos = color_cam_pos / color_cam_pos(3);
			color_cam_pos(3) = 1.0;

			// world space
			Eigen::Vector4f world_pos = color_extr_inv * color_cam_pos;

			world_pos = world_pos / world_pos(3);
			world_pos(3) = 1.0;

			pts_list.push_back(world_pos.head<3>());
			pts_rgb_list.push_back(rgb);
		}
	}

	Eigen::MatrixXf pts(pts_list.size(), 3);
	Eigen::MatrixXf pts_rgb(pts_rgb_list.size(), 3);

	for (size_t i = 0; i < pts_list.size(); i++)
	{

		pts.row(i) = pts_list[i];
		pts_rgb.row(i) = pts_rgb_list[i];

	}

	out_pts = pts;
	out_pts_rgb = pts_rgb;
}

void load_scans(ml::SensorData & sd,
	int max_loaded,
	std::vector<Eigen::MatrixXf> & out_depth,
	std::vector<Eigen::MatrixXf> & out_color,
	std::vector<ml::mat4f> & out_cam_pose,
	std::vector<Eigen::MatrixXf> & out_pcd,
	std::vector<Eigen::MatrixXf> & out_pcd_rgb)
{

	//dimensions of a color/depth frame 
	for (size_t i = 0; i < sd.m_frames.size(); i++)
	{
		//de-compress color and depth values
		unsigned short* depthData = sd.decompressDepthAlloc(i);
		ml::vec3uc* colorData = sd.decompressColorAlloc(i);

		Eigen::MatrixXf depth(sd.m_colorWidth*sd.m_colorHeight, 1);
		Eigen::MatrixXf rgb(sd.m_depthWidth*sd.m_depthHeight, 3);

		// depth
		for (unsigned int i = 0; i < sd.m_depthWidth * sd.m_depthHeight; i++) {
			//convert depth values to m:
			float depth_in_meters = depthData[i] / sd.m_depthShift;

			depth(i) = depth_in_meters;
		}
		out_depth.push_back(depth);

		// color
		for (unsigned int i = 0; i < sd.m_colorWidth * sd.m_colorHeight; i++) {
			//convert depth values to m:
			ml::vec3uc c = colorData[i];
			rgb(i, 0) = c.x / 255.0;
			rgb(i, 1) = c.y / 255.0;
			rgb(i, 2) = c.z / 255.0;
		}
		out_color.push_back(rgb);

		// camera pose
		out_cam_pose.push_back(sd.m_frames[i].getCameraToWorld());


		std::free(colorData);
		std::free(depthData);

		if (i + 1 > max_loaded) {
			break;
		}
	}

	//---------------------------------------------------------
	// unproject

	// camera pose 

	sd.m_calibrationColor.m_intrinsic;
	sd.m_calibrationColor.m_extrinsic;
	//sd.m_calibrationDepth.m_intrinsic;
	//sd.m_calibrationDepth.m_extrinsic;

	Eigen::Matrix4f color_extr;
	Eigen::Matrix4f depth_intr;
	//Eigen::Matrix4f depth_extr;

	ml2eigen(sd.m_calibrationColor.m_extrinsic, color_extr);
	ml2eigen(sd.m_calibrationDepth.m_intrinsic, depth_intr);
	//ml2eigen(sd.m_calibrationDepth.m_extrinsic, depth_extr);

	std::cout << "color_extr=\n" << color_extr << "\n";
	std::cout << "depth_intr=\n" << depth_intr << "\n";

	Eigen::Matrix4f  depth_intr_inv = depth_intr.inverse();
	Eigen::Matrix4f  color_extr_inv = color_extr.inverse();

	out_pcd.clear();
	out_pcd_rgb.clear();

	for (size_t i = 0; i < out_depth.size(); i++)
	{
		Eigen::MatrixXf pts;
		Eigen::MatrixXf pts_rgb;

		unproj(sd.m_depthWidth, sd.m_depthHeight, depth_intr_inv, color_extr_inv, out_depth[i], out_color[i], pts, pts_rgb);

		out_pcd.push_back(pts);
		out_pcd_rgb.push_back(pts_rgb);
	}

	std::cout << "loaded x" << out_pcd.size() << " frames\n";

	//-------------------------------------------------


}

//=================================================================================

int main(int argc, char *argv[])
{	

	g_myctx.mode = 1;
	g_myctx.frame_idx = 0;
	g_myctx.nv_len = 0.5;
	g_myctx.point_size = 2;
	g_myctx.line_width = 2;


	//------------------------------------------------------
	// load scan data
	std::string sens_path = "../data/static_dock.sens";
	if (!igl::file_exists(sens_path))
	{
		std::cout << "[error] scan data not found\nPress any key to exit\n";
		char c;
		std::cin >> c;
		return 1;
	}
	ml::SensorData sd(sens_path);

	load_scans(sd, 3, g_myctx.depths,g_myctx.colors,g_myctx.cam_poses, g_myctx.point_clouds, g_myctx.point_clouds_rgb 	);

	//------------------------------------------------------
	// You task is to finish the smooth_depth function 
	// 

	Eigen::MatrixXf averaged_depth;
	Eigen::MatrixXf averaged_pts;
	Eigen::MatrixXf averaged_pts_rgb;

	average_depth( g_myctx.depths, 100, g_myctx.averaged_depth);

	//------------------------------------------------------
	Eigen::Matrix4f color_extr;
	Eigen::Matrix4f depth_intr;
	ml2eigen(sd.m_calibrationColor.m_extrinsic, color_extr);
	ml2eigen(sd.m_calibrationDepth.m_intrinsic, depth_intr);
	Eigen::Matrix4f  depth_intr_inv = depth_intr.inverse();
	Eigen::Matrix4f  color_extr_inv = color_extr.inverse();

	unproj(sd.m_depthWidth, 
		   sd.m_depthHeight, 
		   depth_intr_inv,  
		   color_extr_inv, 
		   g_myctx.averaged_depth,
		   g_myctx.colors[0], 
		   g_myctx.averaged_pcd,
		   g_myctx.averaged_pcd_rgb);

	//------------------------------------------------------
	// Init the viewer
	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	// Customize the menu
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
		ImGui::SetNextWindowSize(ImVec2(300, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
			"MyProperties", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);

		// point size
		if (ImGui::InputFloat("point_size", &g_myctx.point_size))
		{
			viewer.data().point_size = g_myctx.point_size;
		}

		// frame_index
		if (ImGui::SliderInt("mode", &g_myctx.mode, 1, 2))
		{
			g_myctx.change_mode(viewer);
		}
		ImGui::SameLine(); 
		if (g_myctx.mode == 1) {

			ImGui::Text("1:input, frame=%d",g_myctx.frame_idx);
		}
		else {

			ImGui::Text("2:averaged");
		}

		// frame_index
		if (ImGui::SliderInt("frame_index", &g_myctx.frame_idx, 0, g_myctx.point_clouds.size()))
		{
			g_myctx.reset_display(viewer);
		}
		ImGui::SameLine();  ShowHelpMarker("CTRL+click to input value.");

		// save out pcd.ply
		if (ImGui::Button("save out pcd.ply", ImVec2(-1, 0)))
		{
			std::cout << "save at ./pcd.ply\n";
			Eigen::MatrixXd F;
			igl::writePLY("pcd.ply", g_myctx.point_clouds[g_myctx.frame_idx].cast<double>(),F); 
		}
		
		ImGui::End();
	};


	viewer.callback_key_pressed = &key_press;


	g_myctx.reset_display(viewer);
	g_myctx.set_camera(viewer);
	//------------------------------------------------------
	// Call GUI
	viewer.launch();

}
