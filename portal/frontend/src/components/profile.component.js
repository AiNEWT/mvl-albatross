import React, { Component } from "react";
import { Navigate } from "react-router-dom";
import AuthService from "../services/auth.service";
import UserService from "../services/user.service";

export default class Profile extends Component {
  constructor(props) {
    super(props);

    this.state = {
      redirect: null,
      userReady: false,
      currentUser: { username: "" },
      frame_url: '',
      isServiceRunning: false,
      isLoading: false,
    };
  }

  componentDidMount() {
    const currentUser = AuthService.getCurrentUser();
    if (!currentUser) this.setState({ redirect: "/login" });

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

    this.setState({ currentUser: currentUser, userReady: true })
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

  render() {
    if (this.state.redirect) {
      return <Navigate to={this.state.redirect} />
    }

    const { currentUser, isServiceRunning } = this.state;

    return (
      <div className="container">
        {(this.state.userReady) ?
        <div>
        <header className="jumbotron">
          <h3>
            Launch Albatross Container for {" "}
            <strong>{currentUser.username}</strong>
          </h3>
          <p>
          <strong>Email:</strong>{" "}
          {currentUser.email}
        </p>
        </header>
        </div>
      : null}
      {isServiceRunning ? (
        <div>
          <p>
          Since the container has already been created, you can use it as your home to start using the service.{" "}
          <button onClick={this.createService} style={{ marginTop: '20px' }}>
            Stop Service
          </button>
          </p>
        </div>
          ) : (
        <div>
          <p>
          Container creation can take as long as 1-5 minutes. After 5 minutes, click the Home button to use the service. {" "}
          <button onClick={this.createService} style={{ marginTop: '20px' }}>
            Create Service
          </button>
          </p>
          </div>
          )
      }
    </div>
    );
  }
}
