/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <map>

#include "particle_filter.h"


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 5;
	//particles.resize(num_particles);

	default_random_engine gen;

	double std_x     = std[0];
	double std_y     = std[1];
	double std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for (int i = 0; i < num_particles; i++)
	{
		
		Particle new_particle;
		new_particle.id = i;
		new_particle.x = dist_x(gen);
		new_particle.y = dist_y(gen);
		new_particle.theta = dist_theta(gen);
		new_particle.weight = 1.0;

		particles.push_back(new_particle);

		// cout << "Sample " << i + 1 << " " 
		// 	<< new_particle.x << " " 
		// 	<< new_particle.y << " " 
		// 	<< new_particle.theta << " " 
		// 	<< new_particle.weight << endl;
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	for (int i = 0; i < num_particles; i++)
	{

		// Particle current_particle;
		// current_particle = particles[i];

		Particle & current_particle = particles[i];

		double x_f, y_f, theta_f;
		double x_0, y_0, theta_0;

		x_0     = current_particle.x;
		y_0     = current_particle.y;
		theta_0 = current_particle.theta;

		if (fabs(yaw_rate) > 0.0001)
		{
			x_f     = x_0 + (velocity / yaw_rate) * (sin(theta_0 + (yaw_rate * delta_t)) - sin(theta_0));
			y_f     = y_0 + (velocity / yaw_rate) * (cos(theta_0) - cos(theta_0 + (yaw_rate * delta_t)));
			theta_f = theta_0 + (yaw_rate * delta_t);
		}
		else
		{
			x_f     = x_0 + velocity * delta_t * cos(theta_0);
			y_f     = y_0 + velocity * delta_t * sin(theta_0);
			theta_f = theta_0;
		}

		default_random_engine gen;

		double std_x     = std_pos[0];
		double std_y     = std_pos[1];

		normal_distribution<double> dist_x(x_f, std_x);
		normal_distribution<double> dist_y(y_f, std_y);
		
		// particles[i].x = dist_x(gen);
		// particles[i].y = dist_y(gen);
		// particles[i].theta = theta_f;

		current_particle.x = dist_x(gen);
		current_particle.y = dist_y(gen);
		current_particle.theta = theta_f;

		// cout << "Sample " << i + 1 << " " 
		// 	<< current_particle.x << " " 
		// 	<< current_particle.y << " " 
		// 	<< current_particle.theta << " " 
		// 	<< current_particle.weight << endl;

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	// cout << "std_landmark: " << std_landmark[0] << " " << std_landmark[1] << endl;

	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	double sum_weights = 0.0;

	for (int n = 0; n < num_particles; n++)
	{

		Particle & current_particle = particles[n];

		double x_f, y_f, theta_f, weight;
		
		x_f     = current_particle.x;
		y_f     = current_particle.y;
		theta_f = current_particle.theta;
		// weight  = current_particle.weight;
		weight  = 1.0;


		for (int i = 0; i < observations.size(); i++)
		{
			LandmarkObs obs = observations[i];
			double obs_x, obs_y;

			obs_x = x_f + (obs.x * cos(theta_f)) - (obs.y * sin(theta_f));
			obs_y = y_f + (obs.x * sin(theta_f)) + (obs.y * cos(theta_f));

			// cout << "Particle: " << n << " w: " << weight << " Landmark: " << obs_x << " " << obs_y ;//<< endl;

			vector<double> ranges;

			for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
			{
				Map::single_landmark_s landmark = map_landmarks.landmark_list[j];
				double map_x = landmark.x_f;
				double map_y = landmark.y_f;

				double x_range = map_x - obs_x;
				double y_range = map_y - obs_y;

				// double range   = sqrt(x_range*x_range + y_range*y_range);
				double range   = ((x_range*x_range)/std_x) + ((y_range*y_range)/std_y);

				ranges.push_back(range);
			}

			sort(ranges.begin(), ranges.end());
			// cout << "  ranges: " << ranges[0] << " " << ranges[1] << endl;
			double min_range = ranges[0];

			weight *= exp(-0.5*min_range)/(2*M_PI*std_x*std_y); // seriously??

		}

		weight += 1e-300;
		current_particle.weight = weight;
		sum_weights += weight;
	}

	for (int n = 0; n < num_particles; n++)
	{

		Particle & current_particle = particles[n];

		double weight = current_particle.weight;

		current_particle.weight = weight / sum_weights;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<double> particle_weights;

	for (int n = 0; n < num_particles; n++)
	{

		Particle & current_particle = particles[n];

		double weight = current_particle.weight;

		cout << "Particle: " << n << " w: " << weight << endl;

		particle_weights.push_back(weight);
		
	}


	
	for (int n = 0; n < num_particles; n++)
	{

	}

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
