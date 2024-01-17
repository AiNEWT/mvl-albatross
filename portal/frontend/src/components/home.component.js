import React, { Component } from "react";
import { Navigate } from "react-router-dom";

import AuthService from "../services/auth.service";
import UserService from "../services/user.service";

export default class Home extends Component {
  
  constructor(props) {
    super(props);

    this.state = {
      frame_url: '',
      isServiceRunning: false,
      isLoading: false,
    };
  }

  checkServiceStatus = async () => {
    try {
      UserService.getPublicContent().then(
        response => {
          this.setState({
            frame_url: response.data.frame_url,
            isServiceRunning: response.data.Status === 'running' ? true : false,
            isLoading: response.data.Status === 'running' ? true : false,
          });
        },
        error => {
          this.setState({ redirect: '/login' });
        }
      );
    } catch (error) {
      console.error('Error checking service status:', error);
    }
  };

  createService = async () => {
    this.setState({ isLoading: true });
    UserService.postCreateContainer().then(
      response => {
          //console.log(response.data.container_id)        
          this.setState({ isLoading: false });
      },
      error => {
        this.setState({ redirect: '/login' });
      }
    );
    this.checkServiceStatus();
  };

  componentDidMount() {
    const currentUser = AuthService.getCurrentUser();
    if (!currentUser) this.setState({ redirect: '/login' });
    this.setState({ currentUser: currentUser, userReady: true });

    if (currentUser) {
      try {
        this.checkServiceStatus();
      } catch (error) {
        console.error('Error checking service status:', error);
      }
    } else {
      this.setState({ redirect: '/login' });
    }
  }

  render() {

    if (this.state.redirect) {
      return <Navigate to={this.state.redirect} />
    }

    const {frame_url, isLoading, isServiceRunning } = this.state;

    return (
      <div style={{ textAlign: 'center' }}>
        {!isServiceRunning && isLoading && <p>Loading...</p>}
        {isServiceRunning ? (
          <iframe 
            src={frame_url}// URL of the content to be loaded
            title="iFrameAlbatross"
            width='100%'  // Width of the iframe
            height='1400px'  // Height of the iframe
            frameBorder='0' // Style of the border, set to 0 for no border
            allowFullScreen // Allow full screen if the content supports it
          >
            {/* Fallback content in case iframes are not supported */}
            Your browser does not support iframes.
          </iframe>
          ) : (
            <button onClick={this.createService} style={{ marginTop: '20px' }}>
              Create Service
            </button>
          )
        }
      </div>
    );
  }
}