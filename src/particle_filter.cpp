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

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;

	weights.resize(num_particles);

	random_device rd;
	default_random_engine gen(rd());

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

		Particle current_particle = particles[i];

		double x_f, y_f, theta_f;
		double x_0, y_0, theta_0;

		x_0     = current_particle.x;
		y_0     = current_particle.y;
		theta_0 = current_particle.theta;

		if (fabs(yaw_rate) > 0.00001)
		{
			x_f     = x_0 + (velocity / yaw_rate) * (sin(theta_0 + (yaw_rate * delta_t)) - sin(theta_0));
			y_f     = y_0 + (velocity / yaw_rate) * (cos(theta_0) - cos(theta_0 + (yaw_rate * delta_t)));
			theta_f = theta_0 + (yaw_rate * delta_t);
			if (theta_f> M_PI) theta_f = remainder(theta_f, (2.*M_PI)) - M_PI;
    		if (theta_f<-M_PI) theta_f = remainder(theta_f, (2.*M_PI)) + M_PI;
		}
		else
		{
			x_f     = x_0 + velocity * delta_t * cos(theta_0);
			y_f     = y_0 + velocity * delta_t * sin(theta_0);
			theta_f = theta_0;
		}

		// https://stackoverflow.com/questions/15461140/stddefault-random-engine-generate-values-between-0-0-and-1-0
		random_device rd;
		default_random_engine gen(rd());

		double std_x     = std_pos[0];
		double std_y     = std_pos[1];
		double std_theta = std_pos[2];

		normal_distribution<double> dist_x(0.0, std_x);
		normal_distribution<double> dist_y(0.0, std_y);
		normal_distribution<double> dist_theta(0.0, std_theta);
		
		current_particle.x = x_f + dist_x(gen);
		current_particle.y = y_f + dist_y(gen);
		current_particle.theta = theta_f + dist_theta(gen);

		particles[i] = current_particle;
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

	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	double sum_weights = 0.0;

	for (int n = 0; n < num_particles; n++)
	{

		Particle current_particle = particles[n];

		double x_f, y_f, theta_f, weight;
		
		x_f     = current_particle.x;
		y_f     = current_particle.y;
		theta_f = current_particle.theta;
		weight  = 1.0;

		for (int i = 0; i < observations.size(); i++)
		{
			LandmarkObs obs = observations[i];
			double obs_x, obs_y;

			obs_x = x_f + (obs.x * cos(theta_f)) - (obs.y * sin(theta_f));
			obs_y = y_f + (obs.x * sin(theta_f)) + (obs.y * cos(theta_f));

			vector<double> ranges;

			for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
			{
				Map::single_landmark_s landmark = map_landmarks.landmark_list[j];
				double map_x = landmark.x_f;
				double map_y = landmark.y_f;

				double x_range = map_x - obs_x;
				double y_range = map_y - obs_y;

				double range    = sqrt(x_range*x_range + y_range*y_range);

				ranges.push_back(range);
			}

			auto min_range = min_element(begin(ranges),end(ranges));
			int ix = distance(begin(ranges), min_range);
	
			Map::single_landmark_s closest_landmark = map_landmarks.landmark_list[ix];

			double map_x = closest_landmark.x_f;
			double map_y = closest_landmark.y_f;

			double x_range = map_x - obs_x;
			double y_range = map_y - obs_y;

			double term   = ((x_range*x_range)/(std_x*std_x)) + ((y_range*y_range)/(std_y*std_y));

			weight *= (exp(-0.5*term))/(2*M_PI*std_x*std_y); 
			

		}

		weight += 1e-300;
		current_particle.weight = weight;
		sum_weights += weight;

		particles[n] = current_particle;

	}

	for (int n = 0; n < num_particles; n++)
	{

		Particle & current_particle = particles[n];

		double weight = current_particle.weight;

		current_particle.weight = weight / sum_weights;

		weights[n] = current_particle.weight;


	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<double> particle_weights = weights;

	vector<Particle> resampled_particles;

	// http://www.cplusplus.com/reference/random/uniform_int_distribution/
	default_random_engine generator;
  	uniform_int_distribution<int> distribution(0,num_particles-1);

	int index = distribution(generator);
	double beta = 0.0;
	// https://stackoverflow.com/questions/10158756/using-stdmax-element-on-a-vectordouble
	double mw = *max_element(begin(particle_weights), end(particle_weights));

	for (int n = 0; n < num_particles; n++)
	{
		// http://www.cplusplus.com/reference/random/uniform_real_distribution/
		uniform_real_distribution<double> distribution(0.0,1.0);
		double rnd = distribution(generator);
		beta += rnd * 2.0 * mw;
		while (beta > particle_weights[index])
		{
			beta -= particle_weights[index];
			index = (index + 1) % num_particles;
		}
		Particle particle_copy = particles[index];

		resampled_particles.push_back(particle_copy);
	}

	particles = resampled_particles;

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
